import tempfile, unittest
from pathlib import Path
from fortran_namelist_schema import source_schema, validate_input
REPO=Path('/Users/andersb/Jobb/rslmto_fable_v4b')
class NamelistSchemaSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.schema=source_schema(REPO)
    def test_real_staged_deck_rejects_unknown_control_use_kspace(self):
        p=Path('/private/tmp/c52-c2-stage/decks/c2_u0_block/input.nml')
        with self.assertRaisesRegex(ValueError, "unknown members=.*use_kspace"):
            validate_input(p,self.schema)
    def test_use_kspace_is_valid_only_in_self_group(self):
        text='&self\nuse_kspace=.false.\n/\n&control\nrecur="block"\n/'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.nml'; p.write_text(text)
            self.assertTrue(validate_input(p,self.schema)['schema_valid'])
    def test_unknown_member_rejected(self):
        text='&control\nrecur="block"\nuse_kspace=.false.\n/'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.nml'; p.write_text(text)
            with self.assertRaises(ValueError): validate_input(p,self.schema)
if __name__=='__main__': unittest.main(verbosity=2)
