#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

import unittest
import collections
import types


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


class TestUtilFuncs(unittest.TestCase):
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
