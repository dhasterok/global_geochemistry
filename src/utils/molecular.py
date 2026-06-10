#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 15:55:16 2023

@author: Shavin Kaluthantri
"""

class MolecularWeightCalculator:
    def __init__(self):
        self.element_data = {
            'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182,
            'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
            'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.3050,
            'Al': 26.9815386, 'Si': 28.0855, 'P': 30.973762, 'S': 32.065,
            'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
            'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961,
            'Mn': 54.938045, 'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934,
            'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
            'As': 74.92160, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798,
            'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 'Zr': 91.224,
            'Nb': 92.90638, 'Mo': 95.94, 'Ru': 101.07, 'Rh': 102.90550,
            'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818,
            'Sn': 118.710, 'Sb': 121.760, 'Te': 127.60, 'I': 126.90447,
            'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547,
            'Ce': 140.116, 'Pr': 140.90765, 'Nd': 144.242, 'Sm': 150.36,
            'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.500,
            'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
            'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84,
            'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084,
            'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2,
            'Bi': 208.98040, 'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891
        }

    def molecular_weight(self, formula):
        if isinstance(formula, list):
            return [self._compute_mass(f) for f in formula]
        else:
            return self._compute_mass(formula)

    def _compute_mass(self, formula):
        wt = 0
        ind = 0
        while ind < len(formula):
            mass, ind = self._unit_mass(formula, ind)
            wt += mass
        return wt

    def _unit_mass(self, formula, lind):
        if formula[lind] == '(':
            # Handling ion
            rind = lind + 1
            while rind < len(formula) and formula[rind] != ')':
                rind += 1
            if rind >= len(formula):
                raise ValueError("Error in chemical formula. Closing parenthesis not found.")
            ion = formula[lind + 1:rind]
            mass = self.molecular_weight(ion)
            lind = rind + 1
        else:
            # Handling element
            rind = lind + 1
            while rind < len(formula) and formula[rind].islower():
                rind += 1
            sym = formula[lind:rind]
            mass = self._element_lookup(sym)
            lind = rind

        # Counting the number of atoms
        n = 0
        while lind < len(formula) and formula[lind].isdigit():
            n = 10 * n + int(formula[lind])
            lind += 1
        n = max(n, 1)
        return n * mass, lind

    def _element_lookup(self, sym):
        return self.element_data.get(sym, None)

# Example usage
calculator = MolecularWeightCalculator()
weight = calculator.molecular_weight("Fe2O3")
# print(f"Molecular weight of H2O: {weight} g/mol")
