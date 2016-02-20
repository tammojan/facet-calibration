from __future__ import print_function
import sys
import unittest
import traceback
import os

# Append the path of the module to the syspath
sys.path.append('..')
from jsm_util import get_config


class TestConfig(unittest.TestCase):
    
    def setUp(self):
        # Sorted 1D array
        self.config1 = {"param1": 1, 
                        "param2": {"s1": 2, 
                                   "s3": 3
                                   },
                        "param3": {"default": 4, 
                                   "s2": 5
                                   }
                        }
        self.config2 = {'model_facet': '/data/models/output_2015-12/facet2/Facet_1.sky'}

    def test_param1(self):
        """
        Normal existing parameter
        """
        self.assertEqual(get_config("param1", self.config1), 1)

    def test_param2(self):
        """
        Exising parameter with source
        """
        self.assertEqual(get_config("param2", self.config1, source="s1"), 2)
        self.assertEqual(get_config("param2", self.config1, source="s3"), 3)
        
    def test_param2_no_source(self):
        """
        Existing parameter with other sources and no default
        """
        with self.assertRaises(KeyError):
            get_config("param2", self.config1, source="s2")
        
    def test_param2_no_source_default(self):
        """
        Existing parameter with other sources and default in function
        """
        self.assertEqual(get_config("param2", self.config1, source="s2", default=0), 0)
     
    def test_param3(self):
        """
        Existing parameter with other sources and default in key
        Use exsisting source
        """
        self.assertEqual(get_config("param3", self.config1, source="s2"), 5)
     
    def test_param3_no_source(self):
        """
        Existing parameter with other sources and default in key
        Use default
        """
        self.assertEqual(get_config("param3", self.config1), 4)
    
    def test_param3_no_source2(self):
        """
        Existing parameter with other sources and default in key
        Use default, non-existent source
        """
        self.assertEqual(get_config("param3", self.config1, source="s1"), 4)

    def test_param3_no_source_default(self):
        """
        Existing parameter with other sources and default in key
        Use default over function default when appropriate
        """
        self.assertEqual(get_config("param3", self.config1, default=0), 4)
        self.assertEqual(get_config("param3", self.config1, source="s1", default=0), 4)
        self.assertEqual(get_config("param3", self.config1, source="s2", default=0), 5)
        
    def test_param4(self):
        """
        Test non existing parameter
        """
        self.assertEqual(get_config("param4", self.config1, default=0), 0)
        with self.assertRaises(KeyError):
            get_config("param4", self.config1, source="s2")
        with self.assertRaises(KeyError):
            get_config("param4", self.config1)
    
    def test_real(self):
        """
        Test real case
        """
        self.assertEqual(get_config("model_facet", self.config2, default=""), 
                         '/data/models/output_2015-12/facet2/Facet_1.sky')    

        
if __name__ == '__main__':
    unittest.main()