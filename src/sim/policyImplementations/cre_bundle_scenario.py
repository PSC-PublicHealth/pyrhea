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

import logging

import pyrheautils
from policybase import ScenarioPolicy as BaseScenarioPolicy
from cre_bundle_treatment import CREBundleTreatmentPolicy
from cre_bundle_diagnostic import CREBundleDiagnosticPolicy

_validator = None
_constants_values = '$(CONSTANTS)/cre_bundle_scenario_constants.yaml'
_constants_schema = 'scenario_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

# These have to be in order of start-up time, because it's late and I'm tired.
# Format of tuples is (abbrev, startDate, endDate) where the dates start counting
# after burn-in.  For consistency with the calibration, we assume 11/28/2011 
# corresponds to post-burnin day number 150, and that date is the trigger date for
# this intervention.
# testLocations = [('THC_4058_L', 250, 830),
#                  ('THC_365_L', 320, 830),
#                  ('THC_6130_L', 383, 830),
#                  ('THC_2544_L', 452, 830)
#                  ]
#testLocations = [('THC_4058_L', 4, 6),
#                 ('ADVO_450_H', 4, 7),
#                 ('ADVO_801_H', 4, 7),
#                 ('ADVO_836_H', 5, 8),
#                 ('THC_365_L', 6, 8),
#                 ]
# testLocations = [
#                  ('ABBI_31_S',0,100),
# ('ADAM_119_S',0,100),
# ('ADVE_120_H',0,100),
# ('ADVE_500_H',0,100),
# ('ADVE_5101_H',0,100),
# ('ADVE_701_H',0,100),
# ('ADVO_1775_H',0,100),
# ('ADVO_17800_H',0,100),
# ('ADVO_2320_H',0,100),
# ('ADVO_3435_L',0,100),
# ('ADVO_3815_H',0,100),
# ('ADVO_4440_H',0,100),
# ('ADVO_450_H',0,100),
# ('ADVO_801_H',0,100),
# ('ADVO_836_H',0,100),
# ('ALDE_1221_S',0,100),
# ('ALDE_1420_S',0,100),
# ('ALDE_1535_S',0,100),
# ('ALDE_1545_S',0,100),
# ('ALDE_16450_S',0,100),
# ('ALDE_201_S',0,100),
# ('ALDE_2021_S',0,100),
# ('ALDE_2520_S',0,100),
# ('ALDE_255_S',0,100),
# ('ALDE_275_S',0,100),
# ('ALDE_4660_S',0,100),
# ('ALDE_504_S',0,100),
# ('ALDE_5050_S',0,100),
# ('ALDE_5831_S',0,100),
# ('ALDE_6120_S',0,100),
# ('ALDE_803_S',0,100),
# ('ALDE_820_V',0,100),
# ('ALDE_BOX_S',0,100),
# ('ALEX_800_H',0,100),
# ('ALSH_2840_S',0,100),
# ('AMBA_4900_S',0,100),
# ('APPL_21020_S',0,100),
# ('ARBO_1512_S',0,100),
# ('ARBO_535_S',0,100),
# ('ASTA_134_S',0,100),
# ('ASTO_6300_S',0,100),
# ('ATRI_1425_S',0,100),
# ('AURO_10400_H',0,100),
# ('AURO_1601_S',0,100),
# ('AVEN_4505_S',0,100),
# ('BALL_9300_V',0,100),
# ('BALM_2055_S',0,100),
# ('BAYS_1100_S',0,100),
# ('BEEC_1201_S',0,100),
# ('BELH_11401_S',0,100),
# ('BERK_6909_S',0,100),
# ('BERK_8200_S',0,100),
# ('BETH_3298_S',0,100),
# ('BETH_5025_H',0,100),
# ('BETH_8425_S',0,100),
# ('BEVE_1703_S',0,100),
# ('BIRC_1426_S',0,100),
# ('BMO_10602_S',0,100),
# ('BOUR_133_S',0,100),
# ('BRAD_650_S',0,100),
# ('BRAV_1800_S',0,100),
# ('BRAV_2355_S',0,100),
# ('BRAV_3401_S',0,100),
# ('BRAV_4101_S',0,100),
# ('BRAV_850_S',0,100),
# ('BREN_3705_S',0,100),
# ('BRIA_6800_S',0,100),
# ('BRID_8100_S',0,100),
# ('BRIG_4538_S',0,100),
# ('BRIT_8700_S',0,100),
# ('BROO_1800_S',0,100),
# ('BROO_3506_S',0,100),
# ('BUCK_2625_S',0,100),
# ('BURN_14500_S',0,100),
# ('BUTT_339_S',0,100),
# ('BUTT_431_S',0,100),
# ('BUTT_720_S',0,100),
# ('CALI_2829_S',0,100),
# ('CARL_725_S',0,100),
# ('CCL_2401_S',0,100),
# ('CENT_1401_S',0,100),
# ('CENT_2450_S',0,100),
# ('CENT_25_H',0,100),
# ('CENT_4747_S',0,100),
# ('CHAT_7050_S',0,100),
# ('CHEV_3400_S',0,100),
# ('CHIC_1415_S',0,100),
# ('CHIC_6685_S',0,100),
# ('CHUR_2000_S',0,100),
# ('CITY_5825_S',0,100),
# ('CLAR_1366_S',0,100),
# ('CLAR_150_S',0,100),
# ('CLAR_700_S',0,100),
# ('CLAR_7433_S',0,100),
# ('CLAR_829_S',0,100),
# ('COLU_901_S',0,100),
# ('COMM_1136_S',0,100),
# ('COMM_4314_S',0,100),
# ('COMM_503_S',0,100),
# ('COMM_5645_H',0,100),
# ('CONC_9401_S',0,100),
# ('CONS_1000_S',0,100),
# ('CONT_5336_S',0,100),
# ('COOK_1901_H',0,100),
# ('COPL_2000_H',0,100),
# ('COUN_1125_S',0,100),
# ('COUN_1635_S',0,100),
# ('COUN_2330_S',0,100),
# ('COUN_2406_S',0,100),
# ('COUN_421_S',0,100),
# ('COUN_9700_S',0,100),
# ('COUR_3601_S',0,100),
# ('COVE_2625_S',0,100),
# ('COVE_700_S',0,100),
# ('CRES_13301_S',0,100),
# ('CRES_14255_S',0,100),
# ('CRYS_335_S',0,100),
# ('DECA_353_S',0,100),
# ('DEER_306_S',0,100),
# ('DEKA_2600_S',0,100),
# ('DEKA_2944_S',0,100),
# ('DELN_300_H',0,100),
# ('DOBS_120_S',0,100),
# ('DOLT_14325_S',0,100),
# ('DUPA_400_S',0,100),
# ('EDWA_801_H',0,100),
# ('ELMB_127_S',0,100),
# ('ELMH_200_H',0,100),
# ('ELMH_200_S',0,100),
# ('ELMW_1017_S',0,100),
# ('ELMW_7733_V',0,100),
# ('EMBA_555_S',0,100),
# ('EMER_6801_S',0,100),
# ('EVAN_1300_S',0,100),
# ('EVAN_350_S',0,100),
# ('EVER_10124_S',0,100),
# ('EXCE_5701_S',0,100),
# ('EXTE_8633_S',0,100),
# ('FAIR_222_S',0,100),
# ('FAIR_5061_S',0,100),
# ('FAL_9630_S',0,100),
# ('FORE_6840_S',0,100),
# ('FORT_4621_S',0,100),
# ('FOUN_1000_S',0,100),
# ('FRAN_1055_S',0,100),
# ('FRAN_1201_H',0,100),
# ('FRAN_1270_S',0,100),
# ('FRAN_1423_H',0,100),
# ('FRAN_20201_H',0,100),
# ('FRAN_24_H',0,100),
# ('FRAN_301_H',0,100),
# ('FRAN_4021_S',0,100),
# ('FRAN_40_S',0,100),
# ('FRAN_5454_H',0,100),
# ('FRAN_555_S',0,100),
# ('FRAN_701_H',0,100),
# ('GENE_3856_S',0,100),
# ('GEOR_3623_S',0,100),
# ('GLEN_1511_S',0,100),
# ('GLEN_19330_S',0,100),
# ('GLEN_22660_V',0,100),
# ('GLEN_2451_V',0,100),
# ('GLEN_270_S',0,100),
# ('GLEN_3901_S',0,100),
# ('GLEN_4340_S',0,100),
# ('GLEN_8333_V',0,100),
# ('GOTT_701_H',0,100),
# ('GRAS_4621_S',0,100),
# ('GREE_220_S',0,100),
# ('GROS_6601_S',0,100),
# ('GROV_701_S',0,100),
# ('GROV_9000_S',0,100),
# ('HANC_101_S',0,100),
# ('HANC_203_S',0,100),
# ('HARM_3919_S',0,100),
# ('HEAL_5801_S',0,100),
# ('HEAR_3100_S',0,100),
# ('HEAT_15600_S',0,100),
# ('HELI_1308_S',0,100),
# ('HEND_1700_S',0,100),
# ('HEND_1900_S',0,100),
# ('HEND_251_S',0,100),
# ('HEND_3175_S',0,100),
# ('HEND_8800_S',0,100),
# ('HERI_355_S',0,100),
# ('HERI_5888_S',0,100),
# ('HICK_9246_S',0,100),
# ('HIGH_50_S',0,100),
# ('HILL_1740_S',0,100),
# ('HILL_777_S',0,100),
# ('HOLY_12220_S',0,100),
# ('HOLY_2701_H',0,100),
# ('IMPE_3300_S',0,100),
# ('INDI_1007_H',0,100),
# ('INTE_4815_S',0,100),
# ('JACK_5130_S',0,100),
# ('JACK_7531_H',0,100),
# ('JASP_1104_H',0,100),
# ('JOHN_2901_S',0,100),
# ('JOHN_5909_S',0,100),
# ('JOLI_2230_S',0,100),
# ('JOSE_1127_H',0,100),
# ('KANK_100_S',0,100),
# ('KENS_3405_S',0,100),
# ('KIND_2300_S',0,100),
# ('KIND_3415_S',0,100),
# ('KIND_8400_S',0,100),
# ('KISH_ONE_H',0,100),
# ('LAKE_14718_S',0,100),
# ('LAKE_263_S',0,100),
# ('LAKE_7200_S',0,100),
# ('LAKE_735_S',0,100),
# ('LAKE_7618_S',0,100),
# ('LEMO_12450_S',0,100),
# ('LEXI_10300_S',0,100),
# ('LEXI_14601_S',0,100),
# ('LEXI_165_S',0,100),
# ('LEXI_2100_S',0,100),
# ('LEXI_420_S',0,100),
# ('LEXI_4735_S',0,100),
# ('LEXI_675_S',0,100),
# ('LEXI_730_S',0,100),
# ('LEXI_815_S',0,100),
# ('LEXI_900_S',0,100),
# ('LIBE_610_S',0,100),
# ('LIFE_1000_S',0,100),
# ('LIFE_250_S',0,100),
# ('LINC_960_S',0,100),
# ('LITT_2325_S',0,100),
# ('LITT_2800_H',0,100),
# ('LITT_80_S',0,100),
# ('LONG_1666_S',0,100),
# ('LORE_645_H',0,100),
# ('LOYO_2160_H',0,100),
# ('LUTH_1601_S',0,100),
# ('LUTH_16220_S',0,100),
# ('LUTH_800_S',0,100),
# ('MAJO_10352_S',0,100),
# ('MAJO_2350_S',0,100),
# ('MAJO_3301_S',0,100),
# ('MAJO_4410_S',0,100),
# ('MAJO_5025_S',0,100),
# ('MAJO_601_S',0,100),
# ('MAJO_7935_S',0,100),
# ('MAJO_8380_S',0,100),
# ('MANO_11860_S',0,100),
# ('MANO_1500_S',0,100),
# ('MANO_180_S',0,100),
# ('MANO_1920_S',0,100),
# ('MANO_200_S',0,100),
# ('MANO_2145_S',0,100),
# ('MANO_2773_S',0,100),
# ('MANO_3300_S',0,100),
# ('MANO_4225_S',0,100),
# ('MANO_432_S',0,100),
# ('MANO_512_S',0,100),
# ('MANO_600_S',0,100),
# ('MANO_6300_S',0,100),
# ('MANO_715_S',0,100),
# ('MANO_7850_S',0,100),
# ('MANO_900_S',0,100),
# ('MANO_9401_S',0,100),
# ('MANO_940_S',0,100),
# ('MAPL_50_S',0,100),
# ('MATH_820_S',0,100),
# ('MAYF_5905_S',0,100),
# ('MCAL_18300_S',0,100),
# ('MEMO_3701_H',0,100),
# ('MEMO_527_S',0,100),
# ('MERC_2525_H',0,100),
# ('MERC_901_H',0,100),
# ('MICH_802_S',0,100),
# ('MID_4920_S',0,100),
# ('MIDW_111_S',0,100),
# ('MIDW_1430_H',0,100),
# ('MIDW_1519_S',0,100),
# ('MIDW_8540_S',0,100),
# ('MOME_500_S',0,100),
# ('MONT_5550_S',0,100),
# ('MORR_1095_S',0,100),
# ('MORR_1223_S',0,100),
# ('MORR_150_H',0,100),
# ('MOUN_2028_H',0,100),
# ('MSMC_12935_H',0,100),
# ('MST_120_S',0,100),
# ('MUNS_901_H',0,100),
# ('NILE_9777_S',0,100),
# ('NORR_7001_S',0,100),
# ('NORT_251_H',0,100),
# ('NORT_2650_H',0,100),
# ('NORT_4201_H',0,100),
# ('NORT_660_H',0,100),
# ('NORT_777_H',0,100),
# ('NORT_800_H',0,100),
# ('NORT_9600_H',0,100),
# ('NORW_1044_H',0,100),
# ('NORW_2833_S',0,100),
# ('NORW_6016_S',0,100),
# ('OAK_2013_S',0,100),
# ('OAK_221_S',0,100),
# ('OAK_625_S',0,100),
# ('OAK_9525_V',0,100),
# ('OAKR_323_S',0,100),
# ('OAKT_1660_S',0,100),
# ('OUR_1201_S',0,100),
# ('PALO_10426_S',0,100),
# ('PALO_12251_H',0,100),
# ('PARK_6125_V',0,100),
# ('PARK_665_S',0,100),
# ('PAVI_2217_S',0,100),
# ('PERS_3900_S',0,100),
# ('PETE_310_S',0,100),
# ('PETE_520_S',0,100),
# ('PETE_6141_S',0,100),
# ('PETE_746_S',0,100),
# ('PETE_902_S',0,100),
# ('PINE_1212_S',0,100),
# ('PINN_2222_S',0,100),
# ('PINN_9301_H',0,100),
# ('PLAZ_3249_S',0,100),
# ('PLUM_24_S',0,100),
# ('PLYM_315_S',0,100),
# ('PORT_85_H',0,100),
# ('PRAI_345_S',0,100),
# ('PRES_1001_S',0,100),
# ('PRES_100_L',0,100),
# ('PRES_1100_S',0,100),
# ('PRES_1127_H',0,100),
# ('PRES_1325_H',0,100),
# ('PRES_1700_S',0,100),
# ('PRES_20_S',0,100),
# ('PRES_210_S',0,100),
# ('PRES_333_H',0,100),
# ('PRES_355_H',0,100),
# ('PRES_400_S',0,100),
# ('PRES_480_S',0,100),
# ('PRES_500_H',0,100),
# ('PRES_611_S',0,100),
# ('PRES_6930_S',0,100),
# ('PRES_7435_H',0,100),
# ('PRES_77_H',0,100),
# ('PRES_8001_S',0,100),
# ('PRES_811_S',0,100),
# ('PRES_901_S',0,100),
# ('PROV_1101_S',0,100),
# ('PROV_13259_S',0,100),
# ('PROV_3450_S',0,100),
# ('REGE_6631_S',0,100),
# ('RENA_10935_S',0,100),
# ('RENA_2940_S',0,100),
# ('RENA_4437_S',0,100),
# ('RENS_1309_S',0,100),
# ('REST_16300_S',0,100),
# ('RESU_2380_S',0,100),
# ('RESU_500_S',0,100),
# ('RIDG_12550_S',0,100),
# ('RIDG_6450_S',0,100),
# ('RIVE_1601_S',0,100),
# ('RIVE_350_H',0,100),
# ('RIVI_490_S',0,100),
# ('RML_5601_L',0,100),
# ('ROSE_45_H',0,100),
# ('RREM_2155_S',0,100),
# ('RUSH_1653_H',0,100),
# ('RUSH_520_H',0,100),
# ('SAIN_2875_H',0,100),
# ('SALE_1314_S',0,100),
# ('SENE_1301_S',0,100),
# ('SH_700_S',0,100),
# ('SHAB_409_S',0,100),
# ('SHER_1425_H',0,100),
# ('SHER_1950_S',0,100),
# ('SHER_2534_S',0,100),
# ('SHER_5838_S',0,100),
# ('SHER_7350_S',0,100),
# ('SILV_1900_H',0,100),
# ('SKOK_9615_S',0,100),
# ('SLL_7000_S',0,100),
# ('SLOV_3615_S',0,100),
# ('SNH_1200_S',0,100),
# ('SOUT_1010_S',0,100),
# ('SOUT_19000_S',0,100),
# ('SOUT_2649_S',0,100),
# ('SOUT_3311_S',0,100),
# ('SOUT_8012_H',0,100),
# ('SSC_10330_S',0,100),
# ('SSC_2901_S',0,100),
# ('ST_1555_H',0,100),
# ('ST_326_H',0,100),
# ('ST_3800_S',0,100),
# ('ST_9244_S',0,100),
# ('STA_1725_V',0,100),
# ('STC_4321_H',0,100),
# ('STM_1500_H',0,100),
# ('STP_1400_S',0,100),
# ('SWED_5145_H',0,100),
# ('TABO_1347_S',0,100),
# ('THC_225_L',0,100),
# ('THC_2544_L',0,100),
# ('THC_365_L',0,100),
# ('THC_4058_L',0,100),
# ('THC_6130_L',0,100),
# ('THE_110_S',0,100),
# ('THE_1615_S',0,100),
# ('THE_1740_H',0,100),
# ('THE_1_H',0,100),
# ('THE_2425_S',0,100),
# ('THE_2732_S',0,100),
# ('THE_4390_S',0,100),
# ('THE_4600_S',0,100),
# ('THE_55_S',0,100),
# ('THE_5841_H',0,100),
# ('THE_6000_S',0,100),
# ('THE_600_H',0,100),
# ('THE_606_S',0,100),
# ('THE_710_S',0,100),
# ('THE_908_S',0,100),
# ('THI_5400_S',0,100),
# ('THOR_160_S',0,100),
# ('THOR_850_H',0,100),
# ('TOWE_759_S',0,100),
# ('TRI_2500_S',0,100),
# ('TRIL_1101_S',0,100),
# ('TRIL_1251_S',0,100),
# ('UNIT_6308_H',0,100),
# ('VALL_1302_H',0,100),
# ('VALP_3405_S',0,100),
# ('VHS_1225_H',0,100),
# ('VHS_3249_H',0,100),
# ('VHS_3_H',0,100),
# ('VHS_4646_H',0,100),
# ('VIBR_9509_L',0,100),
# ('WARR_66_S',0,100),
# ('WARR_6700_S',0,100),
# ('WASH_10501_S',0,100),
# ('WASH_11308_S',0,100),
# ('WATE_7445_S',0,100),
# ('WATE_7750_S',0,100),
# ('WAUC_176_S',0,100),
# ('WAUK_1324_H',0,100),
# ('WAUK_919_S',0,100),
# ('WEAL_150_S',0,100),
# ('WEST_311_S',0,100),
# ('WEST_3200_S',0,100),
# ('WEST_6501_S',0,100),
# ('WEST_928_S',0,100),
# ('WHEA_1325_S',0,100),
# ('WHIT_300_S',0,100),
# ('WILL_515_S',0,100),
# ('WILL_546_S',0,100),
# ('WIND_124_S',0,100),
# ('WIND_16000_S',0,100),
# ('WINF_28W141_S',0,100),
# ('WISC_471_S',0,100),
# ('WITT_1200_S',0,100),
# ('WOOD_2242_S',0,100),
# ('WOOD_7444_S',0,100),
# ('WOOD_920_S',0,100),
# ('WRHC_309_S',0,100),
# ('WYND_2180_S',0,100),
# ('ZIKA_5448_S',0,100),
# 
#                  ]

class CREBundleScenario(BaseScenarioPolicy):
    def __init__(self, name, patch):
        super(CREBundleScenario, self).__init__(name, patch)
        self.logThisString = _constants['stringToLogWhenStarting']
        self.locationImplementationInformation = _constants['locationsImplementingScenario']['facilities']
        
        
    def begin(self, callingAgent, timeNow):
        logger.warn(self.logThisString)
        assert hasattr(self.patch, 'allFacilities'), ('patch %s has no list of facilities!'
                                                      % self.patch.name)
        
        #for abbrev, startDate, endDate in testLocations:
        for facImp in self.locationImplementationInformation:
            abbrev = facImp['abbrev']
            startDate = facImp['times']['startDate']
            endDate = facImp['times']['endDate']
            if timeNow != startDate:
                assert(timeNow < startDate), 'It is too late to start intervention at %s' % abbrev
                timeNow = callingAgent.sleep(startDate - timeNow)
                
            for fac in self.patch.allFacilities:
                if fac.abbrev == abbrev:
                    fac.flushCaches()
                    for ward in fac.getWards():
                        ward.iA.flushCaches()
                    if isinstance(fac.diagnosticPolicy, CREBundleDiagnosticPolicy):
                        fac.diagnosticPolicy.setValue('active', True)
                    else:
                        raise RuntimeError('%s does not have a CREBundleDiagnosticPolicy' % abbrev)
                    for tP in fac.treatmentPolicies:
                        if isinstance(tP, CREBundleTreatmentPolicy):
                            tP.setValue('active', True)
                            for ward in fac.getWards():
                                for patient in ward.getPatientList():
                                    tP.initializePatientTreatment(ward,patient)
                            logger.info('Activated CREBundleScenario at %s' % abbrev)
                            break
                    else:
                        raise RuntimeError('%s does not have a CREBundleTreatmentPolicy' % abbrev)
                    break
            else:
                raise RuntimeError('Failed to find the facility %s' % abbrev)
        #for abbrev, startDate, endDate in testLocations:
        for facImp in self.locationImplementationInformation:
            abbrev = facImp['abbrev']
            startDate = facImp['times']['startDate']
            endDate = facImp['times']['endDate']
            print "{0}: {1} {2}".format(abbrev,startDate,endDate)
            
            if timeNow != endDate:
                assert(timeNow < endDate), 'It is too late to start intervention at %s' % abbrev
                timeNow = callingAgent.sleep(endDate - timeNow)
            for fac in self.patch.allFacilities:
                if fac.abbrev == abbrev:
                    fac.flushCaches()
                    for ward in fac.getWards():
                        ward.iA.flushCaches()
                    if isinstance(fac.diagnosticPolicy, CREBundleDiagnosticPolicy):
                        fac.diagnosticPolicy.setValue('active', True)
                    for tP in fac.treatmentPolicies:
                        if isinstance(tP, CREBundleTreatmentPolicy):
                            tP.setValue('active', False)
                            logger.info('Deactivated CREBundleScenario at %s' % abbrev)
                            break
                    else:
                        raise RuntimeError('%s does not have a CREBundleTreatmentPolicy' % abbrev)
                    break
            else:
                raise RuntimeError('Failed to find the facility %s' % abbrev)

def getPolicyClasses():
    return [CREBundleScenario]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)


