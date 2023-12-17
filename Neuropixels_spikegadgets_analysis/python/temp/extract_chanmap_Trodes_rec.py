"""
Extract the probe configuration from the header of spikegadgets raw.rec file.

The ".rec" file format contains:
  * a first  text part with information in an XML structure
  * a second part for the binary buffer

Author: Tatsumi Yoshida
"""

# Import library
import numpy as np
import pandas as pd
from xml.etree import ElementTree
import matplotlib.pyplot as plt
import re
import os

from probeinterface import Probe, ProbeGroup, write_prb
from probeinterface.plotting import plot_probe, plot_probe_group

def get_probegroup_np(local_path):
    header_size = None
    with open(local_path, mode='rb') as f:
        while True:
            line = f.readline()
            if b"</Configuration>" in line:
                header_size = f.tell()
                break
        if header_size is None:
            ValueError("SpikeGadgets: the xml header does not contain '</Configuration>'")

        f.seek(0)
        header_txt = f.read(header_size).decode('utf8')

    # explore xml header
    root = ElementTree.fromstring(header_txt)
    # gconf = sr = root.find('GlobalConfiguration')
    hconf = root.find('HardwareConfiguration')
    sconf = root.find('SpikeConfiguration')

    # Number of channels for each probe
    sources = hconf.findall("SourceOptions")
    num_probe = len(sources)
    num_chans = np.zeros((num_probe))
    for pr in range(num_probe):
        source = sources[pr]
        source_id = source.attrib['name']
        options = source.findall("CustomOption")
        for option in options:
            if option.attrib['name'] == 'channelsOn':
                probe_map = option.attrib['data']
                probe_map = [int(ch) for ch in probe_map.split()]
                num_chans[pr] = sum(probe_map)
    # print(num_chans)

    # Get channel info
    # chan_map = np.zeros((int(sum(num_chans)),5))
    chan_map = []
    for trode in sconf:
        for schan in trode:
            pr = int(schan.attrib['spikeSortingGroup'])
            id = int(trode.attrib['id']) - 1000*(pr+1)
            hwchan = int(schan.attrib['hwChan'])
            x = int(schan.attrib['coord_ml'])
            y = int(schan.attrib['coord_dv'])
            ap = int(schan.attrib['coord_ap'])

            chan_map.append((x,y,pr,id,hwchan))
    chan_map = np.asarray(chan_map)
    chan_map = pd.DataFrame(chan_map,columns=['x','y','probe','id','hwChan'])
    print('Probe configuration:')
    print(chan_map)

    # Generate NP probes
    probes = []
    for pr in range(num_probe):
        positions = chan_map[chan_map['probe']==pr][['x','y']]
        shank_ids = np.ones(int(num_chans[pr])) * pr
        probe = Probe(ndim=2, si_units='um')
        probe.set_contacts(positions=positions, shapes='circle', shape_params={'radius': 5},shank_ids=shank_ids)
        # probe.set_contact_ids(contact_ids=contact_ids)
        probes.append(probe)
    # print(probes)

    # print(probes[0].device_channel_indices)

    # Generate probe groups
    probegroup = ProbeGroup()
    for probe in probes:
        probegroup.add_probe(probe)
    # print('probegroup.get_channel_count()', probegroup.get_channel_count())

    for pr in range(num_probe):
        contact_ids = chan_map[chan_map['probe']==pr]['id'].to_numpy() + 200 * pr
        # print(contact_ids)
        # print(type(contact_ids))
        probes[pr].set_device_channel_indices(contact_ids)
    # print('Device_channel_indices:')
    # print(probegroup.get_global_device_channel_indices())

    plot_probe_group(probegroup, same_axes=False, with_channel_index=False)
    plt.show()

    path_split = local_path.split('\\')
    path_split = path_split[-2].split('.')
    output_path = os.path.dirname(local_path) + '\\' + path_split[0] + '.prb'
    # output_path = path_split[0] + '.prb'
    print('Probe file is written at' + output_path)

    # save to probe file
    write_prb(file=output_path,probegroup=probegroup)
