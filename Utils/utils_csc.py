"""
Utils scripts for utils functions

Modified from: https://github.com/tbardouille/camcan_CSC_beta/blob/master/utils_csc.py
Original Authors: Lindsey Power, Cedric Allain, Thomas Moreau, Alexandre Gramfort, Timothy Bardouille
Modified by: Lindsey Power, August 2024
"""
# %%
import seaborn as sns
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from codecs import ignore_errors
import scipy.signal as ss
from sklearn import cluster
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import pandas as pd
from pathlib import Path
import pickle
from joblib import Memory, Parallel, delayed

import mne
from mne_bids import BIDSPath, read_raw_bids

from alphacsc import BatchCDL, GreedyCDL
from alphacsc.utils.signal import split_signal
from alphacsc.utils.convolution import construct_X_multi

# Global variables #
EXP_PARAMS = {
    "sfreq": 150.,             # 300
    "atom_duration": 0.5,      # 0.5,
    "n_atoms": 20,             # 25,
    "reg": 0.2,                # 0.2,
    "eps": 1e-5,               # 1e-4,
    "tol_z": 1e-3,             # 1e-2
}

CDL_PARAMS = {
    'n_atoms': EXP_PARAMS['n_atoms'],
    'n_times_atom': int(np.round(EXP_PARAMS["atom_duration"] * EXP_PARAMS['sfreq'])),
    'rank1': True, 'uv_constraint': 'separate',
    'window': True,
    'unbiased_z_hat': True,
    'D_init': 'chunk',
    'lmbd_max': 'scaled',
    'reg': EXP_PARAMS['reg'],
    'n_iter': 100,
    'eps': EXP_PARAMS['eps'],
    'solver_z': 'lgcd',
    'solver_z_kwargs': {'tol': EXP_PARAMS['tol_z'],
                        'max_iter': 1000},
    'solver_d': 'alternate_adaptive',
    'solver_d_kwargs': {'max_iter': 300},
    'sort_atoms': True,
    'verbose': 1,
    'random_state': 0,
    'use_batch_cdl': True,
    'n_splits': 10,
    'n_jobs': 5
}

def run_csc(X, **cdl_params):
    """Run a CSC model on a given signal X.

    Parameters
    ----------
    X : numpy.ndarray
        the data to run the CSC on

    cdl_params : dict
        dictionary of CSC parameters, such as 'n_atoms', 'n_times_atoms', etc.

    Returns
    -------
    cdl_model

    z_hat_

    """
    print('Computing CSC')

    cdl_params = dict(cdl_params)
    n_splits = cdl_params.pop('n_splits', 1)
    use_batch_cdl = cdl_params.pop('use_batch_cdl', False)
    if use_batch_cdl:
        cdl_model = BatchCDL(**cdl_params)
    else:
        cdl_model = GreedyCDL(**cdl_params)

    if n_splits > 1:
        X_splits = split_signal(X, n_splits=n_splits, apply_window=True)
        X = X[None, :]
    else:
        X_splits = X.copy()

    # Fit the model and learn rank1 atoms
    print('Running CSC')
    cdl_model.fit(X_splits)

    z_hat_ = cdl_model.transform(X)
    return cdl_model, z_hat_

def get_subject_dipole(subject_id, cdl_model=None, info=None):
    """Compute the atoms' dipoles for a subject for a pre-computed CDL model.

    Parameters
    ----------
    subject_id : str
        the subject id

    cdl_model : alphacsc.ConvolutionalDictionaryLearning instance

    info : mne.Info instance


    Returns
    -------
    dip : mne.Dipole instance

    """
    epochFif, transFif, bemFif = get_paths(subject_id)
    if (cdl_model is None) or (info is None):
        # get participant CSC results
        file_name = RESULTS_DIR / subject_id / get_cdl_pickle_name()
        if not file_name.exists():
            print(f"No such file or directory: {file_name}")
            return
        # load CSC results
        cdl_model, info, _, _ = pickle.load(open(file_name, "rb"))
    # select only grad channels
    meg_indices = mne.pick_types(info, meg='grad')
    info = mne.pick_info(info, meg_indices)
    # compute noise covariance
    cov = mne.make_ad_hoc_cov(info)
    u_hat_ = cdl_model.u_hat_
    evoked = mne.EvokedArray(u_hat_.T, info)
    # compute dipole fit
    dip = mne.fit_dipole(evoked, cov, str(bemFif), str(transFif), n_jobs=6,
                         verbose=False)[0]

    # in DAL code
    # # Fit a dipole for each atom
    # # Read in epochs for task data
    # epochs = mne.read_epochs(epochFif)
    # epochs.pick_types(meg='grad')
    # cov = mne.compute_covariance(epochs)
    # info = epochs.info

    # # Make an evoked object with all atom topographies for dipole fitting
    # evoked = mne.EvokedArray(cdl_model.u_hat_.T, info)

    # # Fit dipoles
    # dip = mne.fit_dipole(evoked, cov, bemFif, transFif, verbose=False)[0]

    return dip


def flip_v(v):
    """Ensure the temporal pattern v peak is positive for each atom.

      If necessary, multiply both u and v by -1.

    Parameter
    ---------
    v: array (n_atoms, n_times_atom)
        temporal pattern

    Return
    ------
    v: array (n_atoms, n_times_atom)
    """

    index_array = np.argmax(np.absolute(v), axis=1)
    val_index = np.take_along_axis(v, np.expand_dims(
        index_array, axis=-1), axis=-1).squeeze(axis=-1)
    v[val_index < 0] *= -1
    return v


def get_atoms_info(subject_id, results_dir=RESULTS_DIR):
    """For a given subject, return a list of dictionary containing all atoms'
    informations (subject info, u and v vectors, dipole informations, changes
    in activation before and after button press).

    Parameters
    ----------
    subject_id : str
        the subject id

    results_dir : Pathlib instance
        Path to all participants CSC pickled results

    Returns
    -------
    new_rows : list of dict
    """

    # get participant CSC results
    file_name = results_dir / subject_id / get_cdl_pickle_name()
    if not file_name.exists():
        print(f"No such file or directory: {file_name}")
        return

    # load CSC results
    cdl_model, info, allZ, _ = pickle.load(open(file_name, "rb"))

    # get informations about the subject
    age, sex, hand = get_subject_info(subject_id, PARTICIPANTS_FILE)
    base_row = {'subject_id': subject_id, 'age': age, 'sex': sex, 'hand': hand}
    # get informations about atoms
    dip = get_subject_dipole(subject_id, cdl_model, info=info)

    new_rows = []
    for kk, (u, v) in enumerate(zip(cdl_model.u_hat_, flip_v(cdl_model.v_hat_))):
        gof, pos, ori = dip.gof[kk], dip.pos[kk], dip.ori[kk]

        # calculate the percent change in activation between different phases of movement
        # -1.25 to -0.25 sec (150 samples)
        pre_sum = np.sum(allZ[:, kk, 68:218])
        # -0.25 to 0.25 sec (75 samples)
        move_sum = np.sum(allZ[:, kk, 218:293])
        # 0.25 to 1.25 sec (150 samples)
        post_sum = np.sum(allZ[:, kk, 293:443])

        # multiply by 2 for movement phase because there are half as many samples
        z1 = (pre_sum - 2 * move_sum) / pre_sum
        z2 = (post_sum - 2 * move_sum) / post_sum
        z3 = (post_sum - pre_sum) / post_sum

        new_rows.append({
            **base_row, 'atom_id': kk, 'u_hat': u, 'v_hat': v, 'dipole_gof': gof,
            'dipole_pos_x': pos[0], 'dipole_pos_y': pos[1], 'dipole_pos_z': pos[2],
            'dipole_ori_x': ori[0], 'dipole_ori_y': ori[1], 'dipole_ori_z': ori[2],
            'pre-move_change': z1, 'post-move_change': z2, 'post-pre_change': z3,
            'focal': (gof >= 95), 'rebound': (z3 >= 0.1),
            'movement_related': (z1 >= 0. and z2 >= 0.6)
        })

    return new_rows


def get_atom_df(subject_ids=SUBJECT_IDS, results_dir=RESULTS_DIR, save=True):
    """ Create a pandas.DataFrame where each row is an atom, and columns are
    crutial informations, such a the subject id, its u and v vectors as well
    as the participant age and sex.

    Parameters
    ----------
    subject_ids : list of str
        list of subject ids to which we want to collect their atoms' info

    results_dir : Pathlib instance
        Path to all participants CSC pickled results

    save : bool
        if True, save output dataframe as csv
        defaults to True

    Returns
    -------
    pandas.DataFrame
    """

    new_rows = Parallel(n_jobs=N_JOBS, verbose=1)(
        delayed(get_atoms_info)(this_subject_id)
        for this_subject_id in subject_ids)

    df = pd.DataFrame()
    for this_new_row in new_rows:
        df = df.append(this_new_row, ignore_index=True)

    if save:
        df.to_csv(results_dir / 'all_atoms_info.csv')
        pickle.dump(df, open(results_dir / 'all_atoms_info.pkl', "wb"))

    return df


def double_correlation_clustering(atom_df, u_thresh=0.4, v_thresh=0.4,
                                  exclude_subs=None,
                                  output_dir=RESULTS_DIR):
    """

    Parameters
    ----------
    atom_df : pandas DataFrame
        each row is an atom, at least the columns 'subject_id', 'atom_id',
        'u_hat' and 'v_hat'    

    """

    if exclude_subs is not None:
        atom_df = atom_df[~atom_df['subject_id'].isin(
            exclude_subs)].reset_index()

    # Calculate the correlation coefficient between all atoms
    u_coefs = np.corrcoef(np.stack(atom_df['u_hat']))
    v_list = np.stack(atom_df['v_hat'])
    v_coefs = np.reshape([np.max(ss.correlate(v1, v2))
                          for v1 in v_list for v2 in v_list],
                         (v_list.shape[0], v_list.shape[0]))

    group_num = 0

    # Make atom groups array to keep track of the group that each atom belongs to
    atom_groups = pd.DataFrame(
        columns=['subject_id', 'atom_id', 'index', 'group_number'])

    for ii, row in atom_df.iterrows():
        unique = True
        subject_id, atom_id = row.subject_id, row.atom_id
        max_corr, max_group = 0, 0

        # Loops through the existing groups and calculates the atom's average
        # correlation to that group
        for group in range(group_num + 1):
            indx = atom_groups[atom_groups['group_number']
                               == group]['index'].tolist()
            # Find the u vector and correlation coefficient comparing the
            # current atom to each atom in the group
            avg_u = np.mean(abs(np.asarray([u_coefs[ii][jj] for jj in indx])))
            avg_v = np.mean(abs(np.asarray([v_coefs[ii][jj] for jj in indx])))

            # check if this group passes the thresholds
            if (avg_u > u_thresh) & (avg_v > v_thresh):
                unique = False
                # If it does, also check if this is the highest cumulative
                # correlation so far
                if (avg_u + avg_v) > max_corr:
                    max_corr = (avg_u + avg_v)
                    max_group = group

        if unique:
            # If a similar group is not found, a new group is create and the
            # current atom is added to that group
            group_num += 1
            group_dict = {'subject_id': subject_id, 'atom_id': atom_id,
                          'index': ii, 'group_number': group_num}
        else:
            # If the atom was similar to at least one group, sorts it into the
            # group that it had the highest cumulative correlation to
            group_dict = {'subject_id': subject_id, 'atom_id': atom_id,
                          'index': ii, 'group_number': max_group}

        # Add to group dataframe and reset unique boolean
        atom_groups = atom_groups.append(group_dict, ignore_index=True)

    if output_dir is not None:
        # Save atomGroups to dataframe
        csv_dir = output_dir / \
            ('u_' + str(u_thresh) + '_v_' + str(v_thresh) + '_atom_groups.csv')
        atom_groups.to_csv(csv_dir)

    group_summary = atom_groups.groupby('group_number')\
        .agg({'subject_id': 'nunique', 'atom_id': 'count'})\
        .rename(columns={'subject_id': 'nunique_subject_id',
                         'atom_id': 'count_atom_id'})\
        .reset_index()

    return atom_groups, group_summary


def single_subject_exclusion(atom_df, u_thresh=0.8, v_thresh=0.8,
                             n_group_thresh=14, output_dir=RESULTS_DIR):
    """

    """

    def procedure(subject_id):

        atom_groups = double_correlation_clustering(
            atom_df=atom_df[atom_df['subject_id'] == subject_id].reset_index(),
            u_thresh=u_thresh, v_thresh=v_thresh, output_dir=None)

        new_row = {'subject_id': subject_id,
                   'exclude': False,
                   'n_groups': atom_groups['group_number'].nunique()}

        if new_row['n_groups'] < n_group_thresh:
            new_row['exclude'] = True

        return new_row

    new_rows = Parallel(n_jobs=N_JOBS, verbose=1)(
        delayed(procedure)(this_subject_id)
        for this_subject_id in np.unique(atom_df.subject_id))

    df = pd.DataFrame()
    for this_new_row in new_rows:
        df = df.append(this_new_row, ignore_index=True)

    if output_dir is not None:
        # Save atomGroups to dataframe
        df.to_csv(output_dir / 'df_single_subject_exclusion.csv')

    return df


def correlation_clustering_atoms(atom_df, threshold=0.4,
                                 output_dir=RESULTS_DIR):
    """

    Parameters
    ----------
    threshold : float
        threshold to create new groups

    Returns
    -------
    groupSummary, atomGroups (and save them XXX)

    """

    # XXX exclude 'bad' subjects (single slustering operation)

    # XXX make it read a pre-saved file
    exclude_subs = ['CC420061', 'CC121397', 'CC420396', 'CC420348', 'CC320850',
                    'CC410325', 'CC121428', 'CC110182', 'CC420167', 'CC420261',
                    'CC322186', 'CC220610', 'CC221209', 'CC220506', 'CC110037',
                    'CC510043', 'CC621642', 'CC521040', 'CC610052', 'CC520517',
                    'CC610469', 'CC720497', 'CC610292', 'CC620129', 'CC620490']

    atom_df = atom_df[~atom_df['subject_id'].isin(exclude_subs)].reset_index()
    # Calculate the correlation coefficient between all atoms
    # u_vector_list = np.asarray(atom_df['u_hat'].values)
    # v_vector_list = np.asarray(atom_df['v_hat'].values)

    # v_coefs = []
    # for v in v_vector_list:
    #     for v2 in v_vector_list:
    #         coef = np.max(ss.correlate(v, v2))
    #         v_coefs.append(coef)
    # v_coefs = np.asarray(v_coefs)
    # v_coefs = np.reshape(v_coefs, (10760, 10760))

    # u_coefs = np.corrcoef(u_vector_list, u_vector_list)[0:10760][0:10760]

    # Calculate the correlation coefficient between all atoms
    u_coefs = np.corrcoef(np.stack(atom_df['u_hat']))
    v_list = np.stack(atom_df['v_hat'])
    v_coefs = np.reshape([np.max(ss.correlate(v1, v2))
                          for v1 in v_list for v2 in v_list],
                         (v_list.shape[0], v_list.shape[0]))

    threshold_summary = pd.DataFrame(
        columns=['Threshold', 'Number of Groups', 'Number of Top Groups'])

    # Set parameters
    u_thresh = threshold
    v_thresh = threshold

    atomNum = 0
    groupNum = 0
    unique = True

    # Make atom groups array to keep track of the group that each atom belongs to
    atomGroups = pd.DataFrame(
        columns=['subject_id', 'atom_id', 'Index', 'Group number'])

    for ii, row in atom_df.iterrows():
        # print(row)
        subject_id, atom_id = row.subject_id, row.atom_id

        max_corr = 0
        max_group = 0

        # Loops through the existing groups and calculates the atom's average
        # correlation to that group
        for group in range(0, groupNum + 1):
            gr_atoms = atomGroups[atomGroups['Group number'] == group]
            inds = gr_atoms['Index'].tolist()
            u_groups = []
            v_groups = []

            # Find the u vector and correlation coefficient comparing the
            # current atom to each atom in the group
            for ind2 in inds:
                u_coef = u_coefs[ii][ind2]
                u_groups.append(u_coef)

                v_coef = v_coefs[ii][ind2]
                v_groups.append(v_coef)

            # average across u and psd correlation coefficients in that group
            u_groups = abs(np.asarray(u_groups))
            avg_u = np.mean(u_groups)

            v_groups = abs(np.asarray(v_groups))
            avg_v = np.mean(v_groups)

            # check if this group passes the thresholds
            if (avg_u > u_thresh) & (avg_v > v_thresh):
                unique = False
                # If it does, also check if this is the highest cumulative
                # correlation so far
                if (avg_u + avg_v) > max_corr:
                    max_corr = (avg_u + avg_v)
                    max_group = group

        # If the atom was similar to at least one group, sorts it into the
        # group that it had the highest cumulative correlation to
        if (unique == False):
            groupDict = {'subject_id': subject_id, 'atom_id': atom_id,
                         'Index': ii, 'Group number': max_group}

        # If a similar group is not found, a new group is create and the
        # current atom is added to that group
        elif (unique == True):
            groupNum += 1
            print(groupNum)
            groupDict = {'subject_id': subject_id, 'atom_id': atom_id,
                         'Index': ii, 'Group number': groupNum}

        # Add to group dataframe and reset unique boolean
        atomGroups = atomGroups.append(groupDict, ignore_index=True)
        unique = True

        # Summary statistics for the current dataframe:

    # Number of distinct groups
    groups = atomGroups['Group number'].tolist()
    groups = np.asarray(groups)
    numGroups = len(np.unique(groups))

    # Number of atoms and subjects per group
    numAtoms_list = []
    numSubs_list = []

    for un in np.unique(groups):
        numAtoms = len(np.where(groups == un)[0])
        numAtoms_list.append(numAtoms)

        groupRows = atomGroups[atomGroups['Group number'] == un]
        sub_list = np.asarray(groupRows['subject_id'].tolist())
        numSubs = len(np.unique(sub_list))
        numSubs_list.append(numSubs)

    numAtoms_list = np.asarray(numAtoms_list)
    meanAtoms = np.mean(numAtoms_list)
    stdAtoms = np.std(numAtoms_list)

    print("Number of groups:")
    print(numGroups)

    print("Average number of atoms per group:")
    print(str(meanAtoms) + " +/- " + str(stdAtoms))

    groupSummary = pd.DataFrame(
        columns=['Group Number', 'Number of Atoms', 'Number of Subjects'])
    groupSummary['Group Number'] = np.unique(groups)
    groupSummary['Number of Atoms'] = numAtoms_list
    groupSummary['Number of Subjects'] = numSubs_list

    numSubs_list = np.asarray(numSubs_list)
    topGroups = len(np.where(numSubs_list >= 12)[0])
    threshold_dict = {'Threshold': threshold,
                      'Number of Groups': numGroups,
                      'Number of Top Groups': topGroups}
    threshold_summary = threshold_summary.append(
        threshold_dict, ignore_index=True)

    if output_dir is not None:

        # Save group summary dataframe
        csv_dir = output_dir + \
            'u_' + str(u_thresh) + '_v_' + str(v_thresh) + '_groupSummary.csv'
        groupSummary.to_csv(csv_dir)

        # Save atomGroups to dataframe
        csv_dir = output_dir + \
            'u_' + str(u_thresh) + '_v_' + str(v_thresh) + '_atomGroups.csv'
        atomGroups.to_csv(csv_dir)

    return groupSummary, atomGroups


def culstering_cah_kmeans(df, data_columns='all', n_clusters=6):
    """Compute a CAH and k-means clustering

    """
    if data_columns == 'all':
        data = np.array(df)
    else:
        data = np.array(df[data_columns])
    # CAH clustering
    clustering = AgglomerativeClustering(n_clusters=n_clusters,
                                         affinity='euclidean',
                                         linkage='ward')
    clustering.fit(data)
    df['labels_cah'] = clustering.labels_
    # k-means clustering
    kmeans = cluster.KMeans(n_clusters=n_clusters)
    kmeans.fit(data)
    df['labels_kmeans'] = kmeans.labels_

    return df


def compute_distance_matrix(atom_df):
    """Compute the distance matrix, where

    ..maths:
        M[i,j] = 1 - \frac{\sqrt{corr_u[i,j]^2 + corr_v[i,j]^2}}{\sqrt{2}}
    """

    corr_u = np.corrcoef(np.stack(atom_df['u_hat']))

    v_list = np.stack(atom_df['v_hat'])
    corr_v = np.reshape([np.max(ss.correlate(v1, v2))
                         for v1 in v_list for v2 in v_list],
                        (v_list.shape[0], v_list.shape[0]))

    D = (1 - np.sqrt(corr_u**2 + corr_v**2) / np.sqrt(2)).clip(min=0)
    np.fill_diagonal(D, 0)  # enforce the 0 in the diagonal
    # ensure that D is symetric while keeping enough precision
    D = np.round(D, 6)

    return D


def reconstruct_class_signal(df, results_dir):
    """ Reonstruct the signal for all atoms in the given dataframe

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe where each row is an atom, and has at least
        the folowing columns :
            subject_id : the participant id associated with the atom
            atom_id : the atom id

    results_dir : Pathlib instance
        Path to all participants CSC pickled results

    Returns
    -------
    X : array-like
        the reconstructed signal

    n_times_atom : int
        the minimum number of timestamps per atom, accross all atoms in the
        input dataframe
    """

    Z_temp = []
    D_temp = []
    min_n_times_valid = np.inf
    for subject_id in set(df['subject_id'].values):
        file_name = results_dir / subject_id / get_cdl_pickle_name()
        cdl_model, info, _, _ = pickle.load(open(file_name, "rb"))
        atom_idx = df[df['subject_id'] ==
                      subject_id]['atom_id'].values.astype(int)
        Z_temp.append(cdl_model.z_hat_[:, atom_idx, :])
        min_n_times_valid = min(min_n_times_valid, Z_temp[-1].shape[2])
        D_temp.append(cdl_model.D_hat_[atom_idx, :, :])

    # combine z and d vectors
    Z = Z_temp[0][:, :, :min_n_times_valid]
    D = D_temp[0]
    for this_z, this_d in zip(Z_temp[1:], D_temp[1:]):
        this_z = this_z[:, :, :min_n_times_valid]
        Z = np.concatenate((Z, this_z), axis=1)
        D = np.concatenate((D, this_d), axis=0)

    n_times_atom = D.shape[-1]

    X = construct_X_multi(Z, D)

    return X, n_times_atom


def get_df_mean(df, col_label='Group number', cdl_params=CDL_PARAMS,
                results_dir=RESULTS_DIR, n_jobs=N_JOBS):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        the clustering dataframe where each row is an atom, and has at least
        the folowing columns :
            subject_id : the participant id associated with the atom
            u_hat : the topomap vector of the atom
            v_hat : the temporal pattern of the atom
            col_label : the cluster result

    col_label : str
        the name of the column that contains the cultering result

    cdl_params : dict
        the CDL parameters to use to compute the mean atom.
        By default, use GreedyCDL, to use BatchCDL, ensure that
        cdl_params['use_batch_cdl'] = True

    results_dir : Pathlib instance
        Path to all participants CSC pickled results

    n_jobs : int
        number of concurrently running jobs

    Returns
    -------
    df_mean : pandas.DataFrame
        columns:
            col_label : clustering label
            u_hat, v_hat : spatial and temporal pattern of the mean atom
            z_hat : activation vector of the mean atom
    """

    # ensure that only one recurring pattern will be extracted
    cdl_params.update(n_atoms=1, n_splits=1)

    def procedure(label):
        # Reconstruct signal for a given class
        X, n_times_atom = reconstruct_class_signal(
            df=df[df[col_label] == label], results_dir=results_dir)
        cdl_params['n_times_atom'] = n_times_atom
        cdl_model, z_hat = run_csc(X, **cdl_params)
        # append dataframe
        new_row = {col_label: label,
                   'u_hat': cdl_model.u_hat_[0],
                   'v_hat': cdl_model.v_hat_[0],
                   'z_hat': z_hat,
                   'n_times_atom': n_times_atom}

        return new_row

    new_rows = Parallel(n_jobs=min(n_jobs, df[col_label].nunique()), verbose=1)(
        delayed(procedure)(this_label) for this_label in df[col_label].unique())

    df_mean = pd.DataFrame()
    for new_row in new_rows:
        df_mean = df_mean.append(new_row, ignore_index=True)

    df_mean.rename(columns={col_label: 'label'}, inplace=True)
    df_mean.to_csv(results_dir / 'df_mean_atom.csv')

    return df_mean


def complete_existing_df(atomData, results_dir=RESULTS_DIR):
    """

    """
    atomData = pd.read_csv('atomData.csv')

    atomData.rename(columns={'Subject ID': 'subject_id',
                             'Atom number': 'atom_id'},
                    inplace=True)

    # participants = pd.read_csv("participants.tsv", sep='\t', header=0)
    participants = pd.read_csv(PARTICIPANTS_FILE, sep='\t', header=0)
    participants['subject_id'] = participants['participant_id'].apply(
        lambda x: x[4:])

    columns = ['subject_id', 'atom_id', 'Dipole GOF',
               'Dipole Pos x', 'Dipole Pos y', 'Dipole Pos z',
               'Dipole Ori x', 'Dipole Ori y', 'Dipole Ori z', 'Focal',
               'Pre-Move Change', 'Post-Move Change',
               'Post-Pre Change', 'Movement-related', 'Rebound']

    atom_df_temp = pd.merge(atomData[columns], participants[[
        'subject_id', 'age', 'sex', 'hand']], how="left", on="subject_id")

    subject_dirs = [f for f in results_dir.iterdir() if not f.is_file()]

    df = pd.DataFrame()
    for subject_dir in subject_dirs:
        subject_id = subject_dir.name
        base_row = {'subject_id': subject_id}
        # get participant CSC results
        file_name = results_dir / subject_id / get_cdl_pickle_name()
        if not file_name.exists():
            print(f"No such file or directory: {file_name}")
            continue

        # load CSC results
        cdl_model, _, _, _ = pickle.load(open(file_name, "rb"))

        for kk, (u, v) in enumerate(zip(cdl_model.u_hat_, flip_v(cdl_model.v_hat_))):
            new_row = {**base_row, 'atom_id': int(kk), 'u_hat': u, 'v_hat': v}
            df = df.append(new_row, ignore_index=True)

    atom_df = pd.merge(atom_df_temp, df, how="left",
                       on=["subject_id", "atom_id"])
    atom_df.rename(columns={col: col.lower().replace(' ', '_')
                            for col in atom_df.columns},
                   inplace=True)
    atom_df.to_csv('all_atoms_info.csv')
    pickle.dump(atom_df, open('all_atoms_info.pkl', "wb"))

    return atom_df


# %%

if __name__ == '__main__':
    atomData = pd.read_csv(RESULTS_DIR / 'atomData.csv')
    col_to_drop = [col for col in atomData.columns if 'Unnamed' in col]
    atomData = atomData.drop(columns=col_to_drop)
    all_atoms_info = complete_existing_df(atomData)
    atom_df = all_atoms_info.copy()

    BEM_DIR = DATA_DIR / "camcan-mne/freesurfer"
    TRANS_DIR = DATA_DIR / "camcan-mne/trans"
    TRANS_HALIFAX_DIR == DATA_DIR / "camcan-mne/trans-halifax"

    BEM_FILES = [f.name for f in BEM_DIR.iterdir()]
    TRANS_FILES = [f.name for f in TRANS_DIR.iterdir()]
    TRANS_HALIFAX_FILES = [f.name for f in TRANS_HALIFAX_DIR.iterdir()]
    df_trans = pd.DataFrame()
    for subject_id in atom_df['subject_id'].unique():
        epochFif, transFif, bemFif = get_paths(subject_id)
        new_row = {'subject_id': subject_id,
                   'bem_in_base': bemFif.name.split('-')[0] in BEM_FILES,
                   'trans_in_base': transFif.name in TRANS_FILES,
                   'trans_in_halifax': transFif.name in TRANS_HALIFAX_FILES}
        df_trans = df_trans.append(new_row, ignore_index=True)

    # %%

    # all_atoms_info = pd.read_csv('./all_atoms_info.csv')
    all_atoms_info = pickle.load(open('all_atoms_info.pkl', "rb"))

    exclude_subs = ['CC420061', 'CC121397', 'CC420396', 'CC420348', 'CC320850',
                    'CC410325', 'CC121428', 'CC110182', 'CC420167', 'CC420261',
                    'CC322186', 'CC220610', 'CC221209', 'CC220506', 'CC110037',
                    'CC510043', 'CC621642', 'CC521040', 'CC610052', 'CC520517',
                    'CC610469', 'CC720497', 'CC610292', 'CC620129', 'CC620490']

    atom_groups, group_summuray = double_correlation_clustering(
        all_atoms_info, u_thresh=0.4, v_thresh=0.4, exclude_subs=exclude_subs,
        output_dir=None)

    # %%
    # subject_id = list(set(all_atoms_info['subject_id'].values))[0]
    subject_id = 'CC110037'
    data_cols = ['u_hat', 'v_hat']
    sub_df = all_atoms_info[all_atoms_info['subject_id'] == subject_id]
    n_sensors = len(all_atoms_info['u_hat'].values[0])
    n_times_atom = len(all_atoms_info['v_hat'].values[0])
    X = pd.DataFrame()
    X[[f'u_{i}' for i in range(n_sensors)]] = pd.DataFrame(
        sub_df.u_hat.tolist(), index=sub_df.index)
    # X[[f'v_{i}' for i in range(n_times_atom)]] = pd.DataFrame(
    #     sub_df.v_hat.tolist(), index=sub_df.index)

    # %%

    neigh = NearestNeighbors(n_neighbors=2, metric='correlation')
    nbrs = neigh.fit(X)
    distances, indices = nbrs.kneighbors(X)
    distances = np.sort(distances, axis=0)
    distances = distances[:, 1]
    plt.plot(distances)
    p = 90
    q = int(X.shape[0] * p / 100)
    plt.vlines(q, 0, 1, linestyles='--', label=f'{p}%')
    plt.legend()
    plt.show()

    eps = round(distances[q], 2)
    print(
        f"epsilon choice so that 90% of individuals have their nearest neightbour in less than epsilon: {eps}")

    # %%
    atom_df = all_atoms_info.copy()
    D = compute_distance_matrix(atom_df)
    print(D.min(), D.max())

    # %%
    list_esp = np.linspace(0.1, 0.5, 9)
    list_min_samples = [1, 2, 3]

    def procedure(eps, min_samples):

        y_pred = DBSCAN(eps=eps, min_samples=min_samples,
                        metric='precomputed').fit_predict(D)
        row = {'eps': eps, 'min_samples': min_samples,
               'n_groups': len(np.unique(y_pred))}

        return row

    new_rows = Parallel(n_jobs=N_JOBS, verbose=1)(
        delayed(procedure)(eps, min_samples)
        for eps in list_esp for min_samples in list_min_samples)

    df_dbscan = pd.DataFrame()
    for new_row in new_rows:
        df_dbscan = df_dbscan.append(new_row, ignore_index=True)

    # %%
    n_groups = df_dbscan.pivot("eps", "min_samples", "n_groups")
    ax = sns.heatmap(n_groups, annot=True)
    ax.set_title('Number of groups obtains with DBScan')
    ax.set_ylabel(r'$\varepsilon$')
    ax.set_xlabel(r"min samples")
    plt.show()
    # %%
