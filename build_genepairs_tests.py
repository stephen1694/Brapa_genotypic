import unittest
import os
from build_genepairs import Ortho
from build_genepairs import Ortho, OrthoDict, parse_args


class Testing(unittest.TestCase):
    def test_ortho(self):
        line = "othro\tAT1.1,AT2\tBR1,BR2,BR3"
        ortho_line = Ortho(line)
        self.assertEqual(ortho_line.subject, ["AT1.1", "AT2"])

    def test_iso(self):
        ortho_dict = OrthoDict()
        line = "othro\tAT1.1,AT2\tBR1,BR2,BR3, BR3"
        ortho_line = Ortho(line)
        ortho_dict.build_dict(ortho_line, drop_dot=False)
        key_set = set(ortho_dict.subject_dict.keys())
        value_set = ortho_dict.subject_dict["AT1.1"]
        self.assertEqual(key_set, {"AT1.1", "AT2"})
        self.assertTrue(value_set, {"BR1, BR2, BR3"})

    def test_noiso(self):
        ortho_dict = OrthoDict()
        line = "othro\tAT1.1,AT2\tBR1,BR2,BR3, BR3"
        ortho_line = Ortho(line)
        ortho_dict.build_dict(ortho_line, drop_dot=True)
        key_set = set(ortho_dict.subject_dict.keys())
        value_set = ortho_dict.subject_dict["AT1"]
        self.assertEqual(key_set, {"AT1", "AT2"})
        self.assertTrue(value_set, {"BR1, BR2, BR3"})

    def test_file(self):
        with open("test_file", "w") as temp_file:
            print("othro\tgroup1\tgroup2\n", file=temp_file)
            print("othro\tAT1.1,AT2\tBR1,BR2,BR3\n", file=temp_file)
            print("othro\tAT5,AT3, AT2\tBR1,BR2,BR3, BR3\n", file=temp_file)
        ortho_dict = OrthoDict()  # make ortholog dictionary object
        ortho_dict.process_file("test_file")  # process each line of input file
        # ortho_dict.print_dict()
        key_set = set(ortho_dict.subject_dict.keys())
        value_set = ortho_dict.subject_dict["AT5"]
        self.assertEqual(key_set, {"AT1", "AT2", "AT3", "AT5"})
        self.assertTrue(value_set, {"BR1, BR2, BR3"})
        os.remove("test_file")

    def test_parser(self):
        args = parse_args(["--keep_isoforms"])
        self.assertTrue(args.keep_isoforms)


if __name__ == "__main__":
    unittest.main()
