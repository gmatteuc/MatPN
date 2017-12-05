# -*- coding: utf-8 -*-
"""
@who: Niccolò Bonacchi (nico)
@when: Created on Thu Dec 17 10:32:29 2015
@where: Mainen Lab - Champalimaud Neuroscience Programme, Lisbon, Portugal
@what: Functions that interact with *.kwik files from Klustaviewa
@why: Because rasters and PSTH's
@how:(code below!)
"""

import pandas as pd


def get_groups(file):
    """
    Extracts all group names and group_ids from *.kwik files.
    Arguments:
        file:   str         (path for *.kwik file)

    Returns:
        dict({group_name: group_id})
            group_:    str
            group_id:      int
    """
    # Open the file
    store = pd.HDFStore(file)
    # Get the cluster groups and names from the cluster_group node
    cluster_group_node = store.get_node('/channel_groups/0/cluster_groups/main')
    cluster_group_mapping = dict((gr._v_attrs.name.decode("utf-8").lower(),
                                 int(gr._v_name)) for gr in cluster_group_node)
    return cluster_group_mapping


def get_time_samples(file, group):
    """
    Gets time samples from Klustaviewa *.kiwk file sorted by cluster identity.

    Arguments:
        file:   str         (path for *.kwik file)
        group:  int or str  (0, 1, 2, 3, ... or
                             'noise', 'MUA', 'good', 'unsorted', custom...)
                             case insensitive
                str    (noise, MUA, good, unsorted, custom_group)

    Returns: pandas.DataFrame object
             [neuron0, neuron1]
           0 [ts0    , ts0]
           1 [ts1    , ts1]
           2 [ts2    , ts2]
           ...
# TODO: check if custom groups are incremental (should work anyway...)
    """
    # Open the file
    store = pd.HDFStore(file)
    # Get the cluster groups and names from the cluster_group node
    cluster_group_mapping = get_groups(file)
    # Check if user is using the group name or its id
    if isinstance(group, str):
        group = cluster_group_mapping[group.lower()]
    elif isinstance(group, int):
        try:
            assert group in cluster_group_mapping.values(), "Don't know the group"
        except IndexError:
            print("Can't find the cluster group id: {}".format(group))
            return
    # Get all the clusters
    clusters = store.get_node('/channel_groups/0/clusters/main')
    # Initialize cluster_group dict with all group_ids found
    cluster_groups = dict.fromkeys(cluster_group_mapping.values())
    # Populate the dict
    for k in cluster_groups.keys():
        cluster_groups.update({k: [int(x._v_name) for x in clusters
                              if x._v_attrs.cluster_group == k]})

    # Get all time_samplees and cluster ids in chronological order
    time_samples = store.get_node('/channel_groups/0/spikes/time_samples').read()
    cluster_ids = store.get_node('/channel_groups/0/spikes/clusters/main').read()

    # Make dataframe of neurons and time_samples
    out = pd.concat([pd.Series(time_samples[cluster_ids == clu], name=clu)
                    for clu in cluster_groups[group]], axis=1)

    return out


if __name__ == '__main__':
    file = '/media/nico/NAS_home/Recording Rig/recordings/pb018/pb018_151120a/t2/datafile001.kwik'
    ts = get_time_samples(file, 2)
    #    cluster_names = [int(x) for x in clusters.__members__]
    #    time_samples_df = pd.DataFrame(time_samples, index=cluster_ids)
