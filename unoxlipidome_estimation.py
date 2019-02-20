# -*- coding: utf-8 -*-
#
# Copyright (C) 2018-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Maria Fedorova
# [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#
# For more info please contact:
#     Developer Zhixu Ni: zhixu.ni@uni-leipzig.de
from itertools import product

import pandas as pd
from scipy.special import comb


class TheoLipidome(object):
    """
    Calculate total number of unoxLipid from given list of fatty acids.
    It can be modified to predict number of all site specific species by using site_specific=True.
    This will use a modified permutation algorithm to make all modifications repeatable at all sites
    while remove mirrored products from TG and cardiolipin.
    LysoPL, MG, and DG also treated with special care for positions
    """

    def __init__(self, fa_lst_path, x_dct):
        """
        Load default settings, see the __main__ function to modify the default values.

        :param fa_lst_path: file path to default FA_list.xlsx file
        :type fa_lst_path: str
        :param x_dct: default values for the number of lipid classes that have n FA residues
        :type x_dct: dict
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
        print('All settings loaded...')

    def get_estimation(self, site_specific=False):
        """
        Calculate Number of sum product for Lipids from FA list

        The number calculated do NOT contain any unmodified lipids.

        :param site_specific: set to False to use combinations only. set to True to generate all site specific species
        :type site_specific: bool
        :return: Number Lipids with all FA chain got oxidation T[all]_ox
        :rtype: int
        """

        tot_lipid = 0
        if site_specific is False:
            tot_lipid += len(self.x_dct['x1']) * self.f
            tot_lipid += len(self.x_dct['x2']) * comb(self.f + 1, 2)
            tot_lipid += len(self.x_dct['x3']) * comb(self.f + 2, 3)
            tot_lipid += len(self.x_dct['x4']) * comb(self.f + 3, 4)
        else:
            tot_lipid += len(self.x_dct['x2']) * (self.f ** 2 - self.f ** 2)
            n_x1 = len(self.x_dct['x1'])
            for x1 in self.x_dct['x1']:
                if x1 in ['LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPS']:
                    n_x1 += 1
                if x1 in ['Monoacylglycerol', 'MG']:
                    n_x1 += 2
            tot_lipid += n_x1 * self.f

            tot_lipid += len(self.x_dct['x2']) * (self.f ** 2)
            # DG with -OH at sn1/sn3 is like PLs and calculated above
            if 'Diacylglycerol' in self.x_dct['x2'] or 'DG' in self.x_dct['x2']:
                tot_lipid += comb(self.f, 2) + self.f
            else:
                pass

            if 'Triacylglycerol' in self.x_dct['x3'] or 'TG' in self.x_dct['x3']:
                tot_lipid += self.get_product_no_mirror(3)
                if len(self.x_dct['x3']) - 1 > 0:
                    tot_lipid += (len(self.x_dct['x3']) - 1) * (self.f ** 3)  # Other Lipid with 3 FA
                else:
                    pass
            else:
                tot_lipid += len(self.x_dct['x3']) * (self.f ** 3 - self.f ** 3)  # Other Lipid with 3 FA

            if 'Cardiolipin' in self.x_dct['x3'] or 'CL' in self.x_dct['x4']:
                tot_lipid += self.get_product_no_mirror(4)
                if len(self.x_dct['x4']) - 1 > 0:
                    tot_lipid += (len(self.x_dct['x4']) - 1) * (self.f ** 4)  # Other Lipid with 4 FA
                else:
                    pass
            else:
                tot_lipid += len(self.x_dct['x4']) * (self.f ** 4)  # Other Lipid with 4 FA

        return int(tot_lipid)

    def get_product_no_mirror(self, n_sn):

        fa_lst = range(self.f)
        all_lst = list(product(fa_lst, repeat=n_sn))
        pre_out_lst = []

        for i in all_lst:
            if tuple(reversed(i)) in pre_out_lst:
                pass
            else:
                pre_out_lst.append(i)

        return len(pre_out_lst)


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

    unoxlipidome = TheoLipidome(usr_fa_lst, usr_lipid_classes)

    #  position non-specific
    print('\nResults for predictions with FA residues combinations only:')
    print('Total number of unoxLipidome: ', unoxlipidome.get_estimation())

    #  position modification type and site specific predictions
    print('\nResults for predictions with sn site specific:')
    print('Total number of unoxLipidome: ', unoxlipidome.get_estimation(site_specific=True))

    print('\nFinished!')
