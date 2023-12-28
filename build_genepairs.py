from typing import Dict, List, Set
import argparse
import sys


class Ortho:
    def __init__(self, line: str):
        _, subject, query = line.strip().split("\t")
        self.subject: List[str] = subject.replace(" ", "").split(",")
        self.query: Set[str] = set(query.replace(" ", "").split(","))


class OrthoDict:
    def __init__(self):
        self.subject_dict: Dict[str, Set[str]] = {}

    def build_dict(self, ortho_line: Ortho, drop_dot: bool = True):
        for subject in ortho_line.subject:
            # remove the isoform id
            if drop_dot:
                subject = subject.split(".")[0]
            # add new set of query genes to subject
            # or join them (union) with existing genes
            if subject in self.subject_dict:
                self.subject_dict[subject] = self.subject_dict[subject].union(
                    ortho_line.query
                )
            else:
                self.subject_dict[subject] = ortho_line.query

    def print_dict(self):
        print(f"subject\tquery")
        for (subject, query) in self.subject_dict.items():
            for q in query:
                print(f"{subject}\t{q}")

    def process_file(self, input_file: str, drop_dot: bool = True):
        with open(input_file) as file:
            _ = file.readline()
            for line in file:
                if line != "\n":
                    ortho_line = Ortho(line)
                    self.build_dict(ortho_line, drop_dot=drop_dot)


def parse_args(args):
    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description="Convert orthofinder to simple one-to-one pairs of gene ids."
    )
    # Required positional argument
    parser.add_argument(
        "orthofinder_file", type=str, help="File from orthofinder to convert"
    )
    # Switch
    parser.add_argument(
        "--keep_isoforms",
        action="store_false",
        default=True,
        help="""Switch to use if the isoform index should be kept. 
        For example, flag will keep '.4' in 'AT234.4', 
        while the default is to remove it.""",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    ortho_dict = OrthoDict()  # make ortholog dictionary object
    ortho_dict.process_file(
        args.orthofinder_file, args.keep_isoforms
    )  # process each line of input file
    ortho_dict.print_dict()  # print contents of ortholog dictionary object
