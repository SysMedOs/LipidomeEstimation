# -*- coding: utf-8 -*-
#
# Copyright (C) 2018-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Maria Fedorova
# [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#
# For more info please contact:
#     Developer Zhixu Ni: zhixu.ni@uni-leipzig.de

import pandas as pd
from scipy.special import comb


class TheoOxLipidome(object):

    """
    Calculate total number of oxLipid from given list of fatty acids.
    Python representation of all equations in Figure 1 in corresponding publication.
    Default values can be modified to calculate user specific unoxlipidome.
    The predicted number refer to the unique combinations of modification types and numbers.
    It can be modified to predict number of all site specific species by using site_specific=True.
    This will use a modified permutation algorithm to make all modifications repeatable at all sites
    while remove mirrored products from TG and cardiolipin.
    LysoPL, MG, and DG also treated with special care for positions
    """

    def __init__(self, fa_lst_path, x_dct, mod_dct):
        """
        Load default settings, see the __main__ function to modify the default values.

        :param fa_lst_path: file path to default FA_list.xlsx file
        :type fa_lst_path: str
        :param x_dct: default values for the number of lipid classes that have n FA residues
        :type x_dct: dict
        :param mod_dct: default values for the modification types
        :type mod_dct: dict
        """
        fa_df = pd.read_excel(fa_lst_path, header=0)
        print('Load FA list:')
        print(fa_df)
        self.f = fa_df.shape[0]  # total number of Free Fatty acids F
        db_lst = fa_df['DB'].values.tolist()
        n_db_lst = list(set(db_lst))
        self.fa_dct = {}
        self.n_lst = []
        for n_db in n_db_lst:
            if n_db > 0:
                self.fa_dct[n_db] = db_lst.count(n_db)
                self.n_lst.append(n_db)
        self.x_dct = x_dct
        self.m_dct = mod_dct
        print('All settings loaded...')

    def get_oap(self, m, site_specific=False):
        """
        Calculate Number of OAP from FA with n_db and m modification sites A\ :sub:`n`\

        The number calculated do NOT contain any unmodified FA

        :param m: number of modification sites
        :type m: int
        :param site_specific: set to False to use combinations only. set to True to generate all site specific species
        :type site_specific: bool
        :return: Number of OAP from FA with n_db and m modification sites A_n
        :rtype: int
        """
        try:
            # return comb((self.m_dct['m_oap'] + 1) + m - 1, m) - 1
            # -1 to remove the unmodified structure
            # simplify the equation above
            if site_specific is False:
                return comb(len(self.m_dct['m_oap']) + m, m) - 1  # A_n
            else:
                return len(self.m_dct['m_oap']) ** m - 1  # -1 to remove the unmodified structure
        except KeyError:
            return False

    def get_ocp(self, m, site_specific=False):
        """
        Calculate Number of OCP from FA with n_db and m modification sites C\ :sub:`n`\

        The number calculated do NOT contain any unmodified FA

        :param m: number of modification sites
        :type m: int
        :param site_specific: set to False to use combinations only. set to True to generate all site specific species
        :type site_specific: bool
        :return: Number of OCP from FA with n_db and m modification sites C_n
        :rtype: int
        """

        m_ocp = len(self.m_dct['m_ocp'])
        if m == 1:
            return m_ocp
        else:
            ocp_count = m_ocp
            m_lst = list(range(2, m))
            if m_lst:
                for i in m_lst:
                    # add back the unmodified segment from oap with n-1 C=C
                    ocp_count += m_ocp * (self.get_oap(i - 1, site_specific=site_specific) + 1)

            return ocp_count  # C_n

    def get_cyclic(self, n_db, site_specific=False):
        """
        Calculate Number of sum product of FA chain got oxidation F_ox
        The number calculated do NOT contain any unmodified FA

        :param n_db: number of C=C bond
        :type n_db: int
        :param site_specific: set to False to use combinations only. set to True to generate all site specific species
        :type site_specific: bool
        :return: Number of Prostane and other products for FA with n_db >= 3 P_n
        :rtype: int
        """
        if n_db >= 3:
            try:
                if site_specific is False:
                    return len(self.m_dct['m_p']) + len(self.m_dct['m_o'])  # P_n
                else:
                    tri_unit_site = n_db - 2  # 3db = 1 unit_site, 5db = 3 unit_sites
                    return tri_unit_site * (len(self.m_dct['m_p']) + len(self.m_dct['m_o']))  # P_n
            except KeyError:
                return 0
        else:
            return 0

    def get_all_oxfa(self, site='bis-allylic', site_specific=False):
        """
        Calculate Number of sum product for FA chain got oxidation F\ :sub:`ox`\

        The number calculated do NOT contain any unmodified FA.

        :param site: modification site mode in {'bis-allylic', 'allylic', 'db'}
        :type site: str
        :param site_specific: set to False to use combinations only. set to True to generate all site specific species
        :type site_specific: bool
        :return: Number oxFA F_ox
        :rtype: int
        """
        tot_fa_ox = 0  # Number of sum product for FA chain got oxidation F_ox

        # Define a shift to calc m from n_db
        m_shift = -1  # set default to use bis-allylic sites only --> m = n_db - 1 --> m_shift = -1
        if site in {'bisallylic', 'bis allylic', 'bis-allylic', 'allylic', 'db', 'n_db', 'C=C', 'n'}:
            if site == 'allylic':
                m_shift = 1
            elif site in {'db', 'n_db', 'C=C', 'n'}:
                m_shift = 0
            else:
                pass  # Use "bis-allylic sites only" mode by default...
        else:
            print('!! number of modification sites not defined !!')
            print('Use "bis-allylic sites only" mode by default...')
            # m_shift = -1

        for n_db in self.n_lst:
            _fa_count = self.fa_dct[n_db]  # Number of fatty acids with n C=C bond B<n>
            # set number of modification sites according to user settings
            m = n_db + m_shift
            _f_ox = _fa_count * (self.get_oap(m, site_specific=site_specific)
                                 + self.get_ocp(m, site_specific=site_specific)
                                 + self.get_cyclic(n_db, site_specific=site_specific)
                                 )
            tot_fa_ox += _f_ox

        return int(tot_fa_ox)

    def get_all_class_1oxfa(self, site='bis-allylic', site_specific=False):
        """
        Calculate Number of sum product for Lipids with max 1 FA chain got oxidation T\ :sub:`ox`\

        The number calculated do NOT contain any unmodified lipids.

        :param site: modification site mode in {'bis-allylic', 'allylic', 'db'}
        :type site: str
        :param site_specific: set to False to use combinations only. set to True to generate all site specific species
        :type site_specific: bool
        :return: Number oxLipids with max 1 FA chain got oxidation T_ox
        :rtype: int
        """
        tot_fa_ox = self.get_all_oxfa(site=site, site_specific=site_specific)  # Number of oxFA F_ox
        tox_one_ox = 0
        if site_specific is False:
            tox_one_ox += len(self.x_dct['x1']) * tot_fa_ox
            tox_one_ox += len(self.x_dct['x2']) * tot_fa_ox * self.f
            tox_one_ox += len(self.x_dct['x3']) * tot_fa_ox * comb(self.f + 1, 2)
            tox_one_ox += len(self.x_dct['x4']) * tot_fa_ox * comb(self.f + 2, 3)
        else:
            n_x1 = len(self.x_dct['x1'])
            for x1 in self.x_dct['x1']:
                if x1 in ['LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPS']:
                    n_x1 += 1
                if x1 in ['Monoacylglycerol', 'MG']:
                    n_x1 += 2
            tox_one_ox += n_x1 * tot_fa_ox

            n_x2 = len(self.x_dct['x2'])
            if 'Diacylglycerol' in self.x_dct['x2'] or 'DG' in self.x_dct['x2']:
                n_x2 += 2  # ox FA can be in sn1/sn3 or sn2
            else:
                pass
            tox_one_ox += n_x2 * tot_fa_ox * self.f

            if 'Triacylglycerol' in self.x_dct['x3'] or 'TG' in self.x_dct['x3']:
                tox_one_ox += (tot_fa_ox * self.f * self.f  # oxFA on sn1/sn3 of TG
                               + tot_fa_ox * (comb(self.f, 2) + self.f)  # oxFA on sn2 of TG
                               )
                if len(self.x_dct['x3']) - 1 > 0:
                    tox_one_ox += 3 * (len(self.x_dct['x3']) - 1) * tot_fa_ox * self.f**2  # Other Lipid with 3 FA
                else:
                    pass
            else:
                tox_one_ox += 3 * len(self.x_dct['x3']) * tot_fa_ox * self.f**2  # Other Lipid with 3 FA

            if 'Cardiolipin' in self.x_dct['x3'] or 'CL' in self.x_dct['x4']:
                tox_one_ox += (tot_fa_ox * self.f**3  # oxFA on sn1/sn4 of CL
                               + tot_fa_ox * self.f**3  # oxFA on sn2/sn3 of TG
                               )
                if len(self.x_dct['x3']) - 1 > 0:
                    tox_one_ox += 4 * (len(self.x_dct['x4']) - 1) * tot_fa_ox * self.f**3  # Other Lipid with 4 FA
                else:
                    pass
            else:
                tox_one_ox += 4 * len(self.x_dct['x4']) * tot_fa_ox * self.f**3  # Other Lipid with 4 FA

        return int(tox_one_ox)  # T_ox

    def get_all_class_alloxfa(self, site='bis-allylic', site_specific=False):
        """
        Calculate Number of sum product for Lipids with ALL FA chain got oxidation T\ :sup:`all`:sub:`ox`\

        The number calculated do NOT contain any unmodified lipids.

        :param site: modification site mode in {'bis-allylic', 'allylic', 'db'}
        :type site: str
        :param site_specific: set to False to use combinations only. set to True to generate all site specific species
        :type site_specific: bool
        :return: Number oxLipids with all FA chain got oxidation T[all]_ox
        :rtype: int
        """
        tot_fa_ox = self.get_all_oxfa(site=site, site_specific=site_specific)  # Number of oxFA F_ox
        tox_all_ox = 0
        if site_specific is False:
            tox_all_ox += len(self.x_dct['x1']) * tot_fa_ox
            tox_all_ox += len(self.x_dct['x2']) * (comb(tot_fa_ox + self.f + 1, 2) - comb(self.f + 1, 2))
            tox_all_ox += len(self.x_dct['x3']) * (comb(tot_fa_ox + self.f + 2, 3) - comb(self.f + 2, 3))
            tox_all_ox += len(self.x_dct['x4']) * (comb(tot_fa_ox + self.f + 3, 4) - comb(self.f + 4, 4))
        else:
            tot_fa_opt = tot_fa_ox + self.f  # all possible FA
            tox_all_ox += len(self.x_dct['x2']) * (tot_fa_opt**2 - self.f**2)
            n_x1 = len(self.x_dct['x1'])
            for x1 in self.x_dct['x1']:
                if x1 in ['LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPS']:
                    n_x1 += 1
                if x1 in ['Monoacylglycerol', 'MG']:
                    n_x1 += 2
            tox_all_ox += n_x1 * tot_fa_ox

            tox_all_ox += len(self.x_dct['x2']) * (tot_fa_opt**2 - self.f**2)
            # DG with -OH at sn1/sn3 is like PLs and calculated above
            if 'Diacylglycerol' in self.x_dct['x2'] or 'DG' in self.x_dct['x2']:
                tox_all_ox += (tot_fa_ox * self.f  # -OH at sn2, 1 oxFA and 1 unox FA at sn1/sn3
                               + tot_fa_ox * tot_fa_ox  # -OH at sn2, 2 oxFA at sn1 + sn3
                               )
            else:
                pass

            if 'Triacylglycerol' in self.x_dct['x3'] or 'TG' in self.x_dct['x3']:
                tox_all_ox += (tot_fa_ox * self.f * self.f  # 1 oxFA on sn1/sn3 of TG
                               + tot_fa_ox * (comb(self.f, 2) + self.f)  # 1 oxFA on sn2 of TG
                               + self.f * (comb(tot_fa_ox, 2) + tot_fa_ox)  # 2 oxFA on sn1 + sn3 of TG
                               + self.f * (tot_fa_ox**2)  # 2 oxFA on sn1 + sn2 / sn2 + sn3 of TG
                               + tot_fa_ox * (comb(tot_fa_ox, 2) + tot_fa_ox)  # 3 oxFA on all sn
                               )
                if len(self.x_dct['x3']) - 1 > 0:
                    tox_all_ox += (len(self.x_dct['x3']) - 1) * (tot_fa_opt**3 - self.f**3)  # Other Lipid with 3 FA
                else:
                    pass
            else:
                tox_all_ox += len(self.x_dct['x3']) * (tot_fa_opt**3 - self.f**3)  # Other Lipid with 3 FA

            if 'Cardiolipin' in self.x_dct['x3'] or 'CL' in self.x_dct['x4']:
                tox_all_ox += (tot_fa_ox * self.f**3  # 1 oxFA on sn1/sn4 of CL
                               + tot_fa_ox * self.f**3  # 1 oxFA on sn2/sn3 of CL
                               # 2 oxFA on sn1 + sn4 of CL
                               + (comb(tot_fa_ox, 2) + tot_fa_ox) * (comb(self.f, 2) + self.f)  
                               # 2 oxFA on sn2 + sn3 of CL
                               + (comb(tot_fa_ox, 2) + tot_fa_ox) * (comb(self.f, 2) + self.f)
                               # 2 oxFA on sn1 + sn3 / sn2 + sn4 of CL
                               + (tot_fa_ox * self.f * tot_fa_ox * self.f)
                               # 2 oxFA on sn1 + sn2 / sn3 + sn4 of CL
                               + (tot_fa_ox * tot_fa_ox * self.f * self.f)
                               # 3 oxFA on sn1 + sn2 + sn3 / sn2 + sn3 + sn4 of CL
                               + (tot_fa_ox * tot_fa_ox * tot_fa_ox * self.f)
                               # 3 oxFA on sn1 + sn3 + sn4 / sn1 + sn2 + sn4 of CL
                               + (tot_fa_ox * self.f * tot_fa_ox * tot_fa_ox)
                               # 4 oxFA on all sn of CL
                               + ((comb(tot_fa_ox, 2) + tot_fa_ox)
                                  * (comb(tot_fa_ox, 2) + tot_fa_ox) - tot_fa_ox)
                               )
                if len(self.x_dct['x4']) - 1 > 0:
                    tox_all_ox += (len(self.x_dct['x4']) - 1) * (tot_fa_opt**4 - self.f**4)  # Other Lipid with 4 FA
                else:
                    pass
            else:
                tox_all_ox += len(self.x_dct['x4']) * (tot_fa_opt**4 - self.f**4)  # Other Lipid with 4 FA

        return int(tox_all_ox)  # T[all]_ox


if __name__ == '__main__':

    # load default FA list
    usr_fa_lst = r'data/FA_list.xlsx'

    usr_lipid_classes = {
        # list of lipid classes with 1 FA
        'x1': ['FA', 'CholesterolEster', 'LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPS',
               'Monoacylglycerol', 'Ceramide', 'Sphingolipid'],
        'x2': ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'Diacylglycerol'],  # list of lipid classes with 2 FA
        'x3': ['Triacylglycerol'],  # list of lipid classes with 3 FA
        'x4': ['Cardiolipin'],  # list of lipid classes with 4 FA
    }

    usr_mod_dct = {
        'm_ocp': ['Aldehyde', 'CarboxylicAcid'],  # list of cleavage terminal
        'm_oap': ['OH', 'OOH', 'KETO', 'EPOXY'],  # list of Oxygen addition modifications
        'm_p': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'],  # list of Prostane Rings
        'm_o': ['D-IsoK', 'E-IsoK', 'TXA', 'TXB'],  # list of Other types
    }

    oxlipidome = TheoOxLipidome(usr_fa_lst, usr_lipid_classes, usr_mod_dct)

    #  position non-specific
    print('\nResults for predictions with oxidation sites = default:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa())
    print('Total number of oxLipid with max 1 oxFA: ', oxlipidome.get_all_class_1oxfa())
    print('Total number of oxLipid with all oxFA: ', oxlipidome.get_all_class_alloxfa())

    print('\nResults for predictions with oxidation sites = number of bis-allylic sites:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa(site='bis-allylic'))
    print('Total number of oxLipid with max 1 oxFA: ', oxlipidome.get_all_class_1oxfa(site='bis-allylic'))
    print('Total number of oxLipid with all oxFA: ', oxlipidome.get_all_class_alloxfa(site='bis-allylic'))

    print('\nResults for predictions with oxidation sites = number of C=C bond:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa(site='db'))
    print('Total number of oxLipid with max 1 oxFA: ', oxlipidome.get_all_class_1oxfa(site='db'))
    print('Total number of oxLipid with all oxFA: ', oxlipidome.get_all_class_alloxfa(site='db'))

    print('\nResults for predictions with oxidation site = number of allylic sites:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa(site='allylic'))
    print('Total number of oxLipid with max 1 oxFA: ', oxlipidome.get_all_class_1oxfa(site='allylic'))
    print('Total number of oxLipid with all oxFA: ', oxlipidome.get_all_class_alloxfa(site='allylic'))

    #  position modification type and site specific predictions
    print('\nResults for predictions with oxidation sites = default + site specific:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa(site_specific=True))
    print('Total number of oxLipid with max 1 oxFA: ', oxlipidome.get_all_class_1oxfa(site_specific=True))
    print('Total number of oxLipid with all oxFA: ', oxlipidome.get_all_class_alloxfa(site_specific=True))

    print('\nResults for predictions with oxidation sites = number of bis-allylic sites + site specific:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa(site='bis-allylic', site_specific=True))
    print('Total number of oxLipid with max 1 oxFA: ',
          oxlipidome.get_all_class_1oxfa(site='bis-allylic', site_specific=True))
    print('Total number of oxLipid with all oxFA: ',
          oxlipidome.get_all_class_alloxfa(site='bis-allylic', site_specific=True))

    print('\nResults for predictions with oxidation sites = number of C=C bond + site specific:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa(site='db', site_specific=True))
    print('Total number of oxLipid with max 1 oxFA: ', oxlipidome.get_all_class_1oxfa(site='db', site_specific=True))
    print('Total number of oxLipid with all oxFA: ', oxlipidome.get_all_class_alloxfa(site='db', site_specific=True))

    print('\nResults for predictions with oxidation site = number of allylic sites + site specific:')
    print('Number of oxFA: ', oxlipidome.get_all_oxfa(site='allylic', site_specific=True))
    print('Total number of oxLipid with max 1 oxFA: ',
          oxlipidome.get_all_class_1oxfa(site='allylic', site_specific=True))
    print('Total number of oxLipid with all oxFA: ',
          oxlipidome.get_all_class_alloxfa(site='allylic', site_specific=True))

    print('\nFinished!')
