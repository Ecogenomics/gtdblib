###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2016'
__credits__ = ['Pierre Chaumeil', 'Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = ['Pierre Chaumeil', 'Donovan Parks']
__email__ = 'uqpchaum@uq.edu.au'

from collections import Counter


"""Trimming Functions used to remove columns with insufficient taxa or poor consensus. ."""


def trim_seqs(seqs, min_per_taxa, consensus, min_per_bp):
    """Trim multiple sequence alignment.

    Adapted from the biolib package.

    Parameters
    ----------
    seqs : d[seq_id] -> sequence
        Aligned sequences.
    min_per_taxa : float
        Minimum percentage of taxa required to retain a column [0,1].
    min_per_bp : float
        Minimum percentage of base pairs required to keep trimmed sequence [0,1].
    Returns
    -------
    dict : d[seq_id] -> sequence
        Dictionary of trimmed sequences.
    dict : d[seq_id] -> sequence
        Dictionary of pruned sequences.
    """

    alignment_length = len(seqs.values()[0])

    # count number of taxa represented in each column
    column_count = [0] * alignment_length
    column_chars = [list() for _ in xrange(alignment_length)]
    for seq in seqs.values():
        for i, ch in enumerate(seq):
            if ch != '.' and ch != '-':
                column_count[i] += 1
                column_chars[i].append(ch)

    mask = [False] * alignment_length
    count_wrong_pa = 0
    count_wrong_cons = 0
    for i, count in enumerate(column_count):
        if count >= min_per_taxa * len(seqs):
            c = Counter(column_chars[i])
            if len(c.most_common(1)) == 0:
                ratio = 0
            else:
                _letter, count = c.most_common(1)[0]
                ratio = float(count) / len(column_chars[i])
            if ratio >= consensus:
                mask[i] = True
            else:
                count_wrong_cons += 1
        else:
            count_wrong_pa += 1

    # trim columns
    output_seqs = {}
    pruned_seqs = {}
    for seq_id, seq in seqs.iteritems():
        masked_seq = ''.join([seq[i] for i in xrange(0, len(mask)) if mask[i]])

        valid_bases = len(masked_seq) - masked_seq.count('.') - masked_seq.count('-')
        if valid_bases < len(masked_seq) * min_per_bp:
            pruned_seqs[seq_id] = masked_seq
            continue

        output_seqs[seq_id] = masked_seq

    return output_seqs, pruned_seqs, count_wrong_pa, count_wrong_cons
