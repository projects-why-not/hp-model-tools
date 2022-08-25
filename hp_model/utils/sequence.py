import Bio.PDB as pdb
import re


class HPSequenceManager:
    __polarity = {
        "ALA": "H",
        "ARG": "P",
        "ASN": "P",
        "ASP": "P",
        "CYS": "P",
        "GLY": "H",
        "GLN": "P",
        "GLU": "P",
        "HIS": "P",
        "ILE": "H",
        "LEU": "H",
        "LYS": "P",
        "MET": "H",
        "PHE": "H",
        "PRO": "H",
        "SER": "P",
        "THR": "P",
        "TRP": "H",
        "TYR": "P",
        "VAL": "H"
    }
    __res_notation = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "",
        "ASP": "",
        "CYS": "",
        "GLY": "",
        "GLN": "",
        "GLU": "",
        "HIS": "",
        "ILE": "",
        "LEU": "",
        "LYS": "",
        "MET": "",
        "PHE": "",
        "PRO": "",
        "SER": "",
        "THR": "",
        "TRP": "",
        "TYR": "",
        "VAL": ""
    }

    def __init__(self, sequence=None, pdb_path=None, pdb_chain=None, hp_sequence=None):
        if sequence is None and pdb_path is None and pdb_chain is None and hp_sequence is None:
            raise Exception("HP Model must be initialized with a protein!")
        if hp_sequence is not None:
            self.__init_with_hp_sequence(hp_sequence)
        elif sequence is not None:
            self.__init_with_sequence(sequence)
        elif pdb_chain is not None:
            self.__init_with_pdb_chain(pdb_chain)
        else:
            self.__init_with_pdb_path(pdb_path)

    def get_hp_sequence(self):
        return self.__hp_sequence

    # INITIALIZATION OPTIONS

    def __init_with_hp_sequence(self, hp_sequence):
        self.__hp_sequence = self.parse_hp_sequence(hp_sequence)

    def __init_with_sequence(self, sequence):
        hp_seq = "".join([self.__polarity[resname] for resname in sequence])
        self.__init_with_hp_sequence(hp_seq)

    def __init_with_pdb_chain(self, chain):
        seq = [res.resname for res in chain.child_list]
        self.__init_with_sequence(seq)

    def __init_with_pdb_path(self, pdb_path):
        # WARN: only appropriate PDBs are expected!

        struct = pdb.PDBParser().get_structure("s", pdb_path)
        chain = struct[0].child_list[0]
        self.__init_with_pdb_chain(chain)

    # SEQUENCE MUTATIONS

    @classmethod
    def parse_hp_sequence(cls, sequence):
        if "^" not in sequence:
            return sequence

        parts = re.findall(r"\(([P,H,\^,0-9*,\(,\)]*)\)\^([0-9]*)|([H,P])(\^[0-9]*)?",
                           sequence)
        print(sequence)
        print("\t", parts)
        ans = "".join([cls.parse_hp_sequence(part[0]) * int(part[1]) if len(part[0]) > 0 else cls.parse_hp_sequence(part[2]) * int(part[3][1:] if len(part[3]) > 0 else 1) for part in parts])

        return ans

        # split into parts: naive



        # print(parts)
