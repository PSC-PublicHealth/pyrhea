#! /usr/bin/env python

###################################################################################
# Copyright   2015, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
# =============================================================================== #
#                                                                                 #
# Permission to use, copy, and modify this software and its documentation without #
# fee for personal use within your organization is hereby granted, provided that  #
# the above copyright notice is preserved in all copies and that the copyright    #
# and this permission notice appear in supporting documentation.  All other       #
# restrictions and obligations are defined in the GNU Affero General Public       #
# License v3 (AGPL-3.0) located at http://www.gnu.org/licenses/agpl-3.0.html  A   #
# copy of the license is also provided in the top level of the source directory,  #
# in the file LICENSE.txt.                                                        #
#                                                                                 #
###################################################################################

_rhea_svn_id_ = "$Id$"

import unittest
import collections
import types
import yaml
import jsonschema
import logging

logger = logging.getLogger(__name__)


def enum(*sequential, **named):
    """
    Thanks to stackoverflow user Alec Thomas for this answer to 'How can I represent an
    'Enum' in Python?'
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['names'] = reverse
    return type('Enum', (), enums)


def namedtuple(typename, field_names, verbose=False, rename=False, field_types=None):
    """
    This is equivalent to collections.namedtuple(), but adds some support for printing
    enum strings and for pickling.
    """
    newTp = collections.namedtuple(typename, field_names, verbose=verbose, rename=rename)

    if field_types is not None:
        assert isinstance(field_types, types.ListType), "field_types is not a list"
        assert len(field_types) == len(field_names), \
            "field_types and field_names do not have matching lengths"

        def newrepr(self):
            'Return a nicely formatted representation string'
            bits = [typename, '(']
            start = True
            for fn, val, tp in zip(field_names, self, field_types):
                if tp is not None and tp.__name__ == 'Enum':
                    val = tp.names[val]
                else:
                    val = repr(val)
                if start:
                    start = False
                    bits.extend([fn, '=', val])
                else:
                    bits.extend([', ', fn, '=', val])
            bits.append(')')
            return ''.join(bits)

        newTp.__repr__ = newrepr

    assert newTp.__name__ not in globals(), \
        "module %s already has a type named %s" % (__name__, newTp.__name__)
    globals()[newTp.__name__] = newTp  # So the pickle can find it later

    return newTp


def importConstants(valuePath, schemaPath):
    with open(valuePath, 'rU') as f:
        cJSON = yaml.safe_load(f)
    with open(schemaPath, 'rU') as f:
        schemaJSON = yaml.safe_load(f)
    validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
    nErrors = sum([1 for e in validator.iter_errors(cJSON)])  # @UnusedVariable
    if nErrors:
        for e in validator.iter_errors(cJSON):
            logger.error('Schema violation: %s: %s' %
                         (' '.join([str(word) for word in e.path]), e.message))
        raise RuntimeError('%s does not satisfy the schema %s' %
                           (valuePath, schemaPath))
    return cJSON


class SingletonMetaClass(type):
    """
    Thanks again to stackoverflow, this is a Singleton metaclass for Python.

    see http://stackoverflow.com/questions/6760685/creating-a-singleton-in-python

    Note that 'self' in this case should properly be 'cls' since it is a class, but
    the syntax checker doesn't like that.
    """
    _instances = {}

    def __call__(self, *args, **kwargs):  # stupid syntax checker does not understand metaclasses
        if self not in self._instances:
            self._instances[self] = super(SingletonMetaClass, self).__call__(*args, **kwargs)
        return self._instances[self]


class TestUtilFuncs(unittest.TestCase):
    def test_singletonmeta(self):
        class MyClass(object):
            __metaclass__ = SingletonMetaClass
            instanceCounter = 0

            def __init__(self):
                self.id = self.instanceCounter
                self.instanceCounter += 1

        inst1 = MyClass()
        inst2 = MyClass()
        self.assertTrue(inst1 is inst2, "Singleton is-ness failed")
        self.assertTrue(inst1.id == inst2.id, "Singleton ids do not match")

    def test_namedtuple1(self):
        CareTier = enum('HOME', 'HOSP')
        PatientStatus = namedtuple('PatientStatus', ['careTier', 'thing', 'age'],
                                   field_types=[CareTier, None, types.IntType])
        thing = PatientStatus(CareTier.HOSP, 'first', 7)
        self.assertTrue(str(thing) == "PatientStatus(careTier=HOSP, thing='first', age=7)",
                        "namedtuple.__str__ failed")
        self.assertTrue(repr(thing) == "PatientStatus(careTier=HOSP, thing='first', age=7)",
                        "namedtuple.__repr__ failed")

    def test_namedtuple2(self):
        import pickle
        CareTier = enum('HOME', 'HOSP')
        NewStatus = namedtuple('NewStatus', ['careTier', 'thing', 'age'],
                               field_types=[CareTier, None, types.IntType])
        thing = NewStatus(CareTier.HOSP, 'first', 7)
        pStr = pickle.dumps(thing)
        thing2 = pickle.loads(pStr)
        self.assertTrue((thing.careTier == thing2.careTier) and (thing.thing == thing2.thing)
                        and (thing.age == thing2.age),
                        "namedtuple pickling failed")
        self.assertTrue(str(thing2) == str(thing), "namedtuple pickling fails to preserve __str__")
