#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 João Pedro Rodrigues
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Fixes insertion codes in a PDB file.
Works by deleting an insertion code and shifting the residue numbering of
downstream residues. Allows for picking specific residues to delete insertion
codes for.
Usage:
    python pdb_fixinsert.py [-<option>] <pdb file>
Example:
    python pdb_fixinsert.py 1CTF.pdb  # delete ALL insertion codes
    python pdb_fixinsert.py -A9,B12 1CTF.pdb  # deletes ins. codes for res
                                              # 9 of chain A and 12 of chain B.
This program is part of the `pdb-tools` suite of utilities and should not be
distributed isolatedly. The `pdb-tools` were created to quickly manipulate PDB
files using the terminal, and can be used sequentially, with one tool streaming
data to another. They are based on old FORTRAN77 code that was taking too much
effort to maintain and compile. RIP.
"""

__author__ = "Joao Rodrigues"
__email__ = "j.p.g.l.m.rodrigues@gmail.com"



def run(fhandle, option_list):
    """
    Delete insertion codes (at specific residues).
    By default, removes ALL insertion codes on ALL residues. Also bumps
    the residue numbering of residues downstream of each insertion.
    This function is a generator.
    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.
    option_list : list
        List of insertion options to act on.
        Example ["A9", "B12"]. An empty list performs the default
        action.
    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    """

    option_set = set(option_list)  # empty if option_list is empty

    # Keep track of residue numbering
    # Keep track of residues read (chain, resname, resid)
    offset = 0
    prev_resi = None
    seen_ids = set()
    clean_icode = False
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
    for line in fhandle:
        if line.startswith(records):
            res_uid = line[17:27]  # resname, chain, resid, icode
            id_res = line[21] + line[22:26].strip()  # A99, B12
            has_icode = line[26].strip()  # ignore ' ' here

            # unfortunately, this is messy but not all PDB files follow a nice
            # order of ' ', 'A', 'B', ... when it comes to insertion codes..
            if prev_resi != res_uid:  # new residue
                # Does it have an insertion code
                # OR have we seen this chain + resid combination before?
                # #2 to catch insertions WITHOUT icode ('A' ... ' ' ... 'B')
                if (has_icode or id_res in seen_ids):
                    # Do something about it
                    # if the user provided options and this residue is in them
                    # OR if the user did not provide options
                    if (not option_set) or (id_res in option_set):
                        clean_icode = True
                    else:
                        clean_icode = False
                else:
                    clean_icode = False

                prev_resi = res_uid

                if id_res in seen_ids:  # offset only if we have seen this res.
                    offset += 1

            if clean_icode:  # remove icode
                line = line[:26] + ' ' + line[27:]

            # Modify resid if necessary
            resid = int(line[22:26]) + offset
            line = line[:22] + str(resid).rjust(4) + line[26:]
            seen_ids.add(id_res)

            # Reset offset on TER
            if line.startswith('TER'):
                offset = 0
        yield line



def fixInsert(pdb_file):
    # Check Input
    pdbfh = open(f"../input/{pdb_file}", 'r')
    option_list = []
    # Do the job
    new_pdb = run(pdbfh, option_list)

    with open(f"../input/fixed_{pdb_file}", "w+") as _file:
        for line in new_pdb:
            _file.write(line)
    pdbfh.close()
    return f"fixed_{pdb_file}"


if __name__ == '__main__':
    fixInsert("6waq.pdb")
