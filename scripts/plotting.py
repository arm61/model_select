"""
Nice plotting.

Copyright (c) Andrew R. McCluskey

Distributed under the terms of the MIT License

@author: Andrew R. McCluskey
"""

# pylint: disable=R0902

import sys
import itertools
import numpy as np
import h5py
import matplotlib.pyplot as plt
from dynesty import plotting, utils
from uncertainties import ufloat
import fig_params

CON_LIST = [
    "d13acmw",
    "d13d2o",
    "hd2o",
    "d70acmw",
    "d70d2o",
    "d83acmw",
    "d83d2o",
]

VARIABLES = [
    "th",
    "mvh",
    "tt", 
    "phit",
    "mvt",
    "rough", 
]


sys.path.append('scripts/')


class Plotting:
    """
    Cleaning up plotting class
    """
    def __init__(self):
        """
        Initialisation.
        """
        self.variables = ["th", "mvh", "tt", "phit", "mvt", "rough"]
        self.latex = [
            r"$d_h$",
            r"$V_h$",
            r"$d_t$",
            r"$\phi_t$",
            r"$V_t$",
            r"$\sigma$",
        ]
        self.units = [
            r"$d_h$/Å",
            r"$V_h$/Å$^3$",
            r"$d_t$/Å",
            r"$\phi_t$",
            r"$V_t$/Å$^3$",
            r"$\sigma$/Å",
        ]
        self.variables_p = []
        self.latex_p = []
        self.units_p = []
        self.logz = None
        self.logz_err = None
        self.best_latex = []
        self.best_ev = []
        self.best_file = None

    def generate_permutations(self):
        """
        Generate all of the permutations of the variables.
        """
        for i in range(len(self.variables)):
            for j in itertools.combinations(self.variables, i + 1):
                self.variables_p.append('_'.join(j))
            for j in itertools.combinations(self.latex, i + 1):
                self.latex_p.append('/'.join(j))
            for j in itertools.combinations(self.units, i + 1):
                self.units_p.append(j)

    def read_all_logz(self):
        """
        Read the evidence and uncertainties from the different files.
        """
        self.logz = np.zeros(len(self.variables_p))
        self.logz_err = np.zeros(len(self.variables_p))
        for i, var in enumerate(self.variables_p):
            try:
                file_h5 = h5py.File('output/model^{}.h5'.format(var), 'r')
                self.logz[i] = file_h5['logz'][-1]
                self.logz_err[i] = file_h5['logzerr'][-1]
            except FileNotFoundError:
                print('File output/model^{}.h5 not found.'.format(var))
        self.best_latex = [
            self.latex_p[np.argsort(self.logz)[-1]],
            self.latex_p[np.argsort(self.logz)[-2]],
        ]
        self.best_ev = [
            ufloat(
                np.sort(self.logz)[-1],
                self.logz_err[np.argsort(self.logz)[-1]],
            ),
            ufloat(
                np.sort(self.logz)[-2],
                self.logz_err[np.argsort(self.logz)[-2]],
            ),
        ]

    def write_latex_logz(self):
        """
        Write evidence and uncertainty to files.
        """
        for i in range(len(self.logz)):
            file_latex = open(
                'paper/{}_ev.txt'.format(self.variables_p[i]),
                'w',
            )
            file_latex.write(
                r'\num{' + r'{:.1f}'.format(
                    self.logz[i],
                ) + r'\pm' + r'{:.1f}'.format(self.logz_err[i]) + r'}',
            )
            file_latex.close()
        file_latex = open('paper/best_label.txt', 'w')
        file_latex.write(self.best_latex[0])
        file_latex.close()

        file_latex = open('paper/next_best_label.txt', 'w')
        file_latex.write(self.best_latex[1])
        file_latex.close()

        file_ev = open('paper/best_ev.txt', 'w')
        file_ev.write(
            r'\num{' + r'{:.1f}'.format(
                self.best_ev[0].n
            ) + r'\pm' + r'{:.1f}'.format(
                self.best_ev[0].s
            ) + r'}',
        )
        file_ev.close()

        file_ev = open('paper/next_best_ev.txt', 'w')
        file_ev.write(
            r'\num{' + r'{:.1f}'.format(
                self.best_ev[1].n
            ) + r'\pm' + r'{:.1f}'.format(
                self.best_ev[1].s
            ) + r'}',
        )
        file_ev.close()

        diff = (self.best_ev[0] - self.best_ev[1]) * 2
        file_diff = open('paper/diff_ev.txt', 'w')
        file_diff.write(
            r'\num{' + r'{:.1f}'.format(
                diff.n,
            ) + r'\pm' + r'{:.1f}'.format(diff.s) + r'}',
        )
        file_diff.close()

    def write_latex_logz_table(self):
        """
        Write the evidence values to a table.
        """
        file_latex = open('paper/ev_table.txt', 'w')
        for i, evidence in enumerate(np.argsort(self.logz)[::-1]):
            file_latex.write(
                r'{}'.format(
                    self.latex_p[evidence],
                ) + r' & \num{' + r'{:.1f}'.format(
                    self.logz[evidence]
                ) + r'\pm' + r'{:.1f}'.format(self.logz_err[evidence]) + r'}',
            )
            if i % 2 == 1:
                file_latex.write(' \\\\ \n')
            else:
                file_latex.write(' & ')
            i += 1
        file_latex.close()

    def plot_logz(self):
        """
        Plot the variation in evidence as a function of the valibles.
        """
        axes = plt.subplots(figsize=(25, 6))[1]
        axes.errorbar(
            range(len(self.logz)),
            self.logz,
            self.logz_err,
            marker='o',
            ls='',
        )
        axes.axhline(self.logz.max(), c='k', ls='--')
        axes.axvline(self.logz.argmax(), c='k', ls='--')
        axes.set_xticks(range(len(self.logz)))
        axes.set_xticklabels(self.latex_p)
        axes.set_xlim((-1, len(self.logz)))
        plt.xticks(rotation=90)
        axes.set_ylim(
            (
                self.logz[self.logz != 0].min()-500,
                self.logz[self.logz != 0].max()+500,
            )
        )
        axes.set_xlabel(r'Free Parameters')
        axes.set_ylabel(r'$\ln\{p(\mathbf{D}|H)\}$')
        plt.savefig('paper/evidence.pdf')
        plt.close()

    def plot_best_per_number(self):
        """
        Plot the best number of free parameters at each number of parameters.
        """
        best_z = np.zeros((len(self.variables)))
        best_z_labels = [''] * 6
        for i in range(len(self.variables)):
            for j in itertools.combinations(self.variables, i + 1):
                if self.logz[self.variables_p.index('_'.join(j))] > best_z[i]:
                    best_z[i] = self.logz[self.variables_p.index('_'.join(j))]
                    best_z_labels[i] = self.latex_p[
                        self.variables_p.index('_'.join(j))
                    ]

        axes = plt.subplots(figsize=(10, 6))[1]
        axes.plot(range(0, 6), best_z, 'o-')
        axes.set_xticks(range(0, 6))
        axes.set_xticklabels(best_z_labels)
        plt.xticks(rotation=90)
        axes.set_xlabel(r'Free Parameters')
        axes.set_ylabel(r'$\ln\{p(\mathbf{D}|H)\}$')
        axes.set_ylim(6.8e3, 6.88e3)
        plt.savefig('paper/best_per.pdf')
        plt.close()

    def plot_correlation_best(self):
        """
        Plot the corner plot for the best evidence.
        """
        self.best_file = self.variables_p[np.argsort(self.logz)[-1]]
        file_h5 = h5py.File('output/model^{}.h5'.format(self.best_file), 'r')
        plt.subplots(figsize=(11.25, 11.25))
        plotting.cornerplot(
            file_h5,
            labels=self.units_p[np.argsort(self.logz)[-1]],
            color=list(fig_params.TABLEAU.values())[0],
        )
        plt.savefig('paper/post_best.pdf')
        plt.close()

    def plot_for_all(self):
        """
        Plot the corner plot where all parameters are varied.
        """
        file_h5 = h5py.File('output/model^th_mvh_tt_phit_mvt_rough.h5', 'r')
        axes = plt.subplots(figsize=(10, 5))[1]
        axes.errorbar(
            range(1, len(file_h5['logz']) + 1),
            file_h5['logz'],
            file_h5['logzerr'],
            c=list(fig_params.TABLEAU.values())[0],
        )
        axes.set_xscale('log')
        axes.set_xlabel(r'Samples')
        axes.set_ylabel(r'$\ln\{p(\mathbf{D}|H)\}$')
        plt.savefig('paper/iterations.pdf')
        plt.close()
        plt.subplots(figsize=(11.25, 11.25))
        plotting.cornerplot(
            file_h5,
            labels=self.units_p[-1],
            color=list(fig_params.TABLEAU.values())[0],
        )
        plt.savefig('paper/post_all.pdf')
        plt.close()

    def plot_for_any(self, vars, i):
        """
        Plot the corner plot where all parameters are varied.
        """
        file_h5 = h5py.File('output/model^{}.h5'.format(vars), 'r')
        plt.subplots(figsize=(11.25, 11.25))
        plotting.cornerplot(
            file_h5,
            labels=self.units_p[i],
            color=list(fig_params.TABLEAU.values())[0],
        )
        plt.savefig('paper/post_{}.pdf'.format(vars))
        plt.close()

    def make_plots_latex(self):
        """
        Make a latex snippet for all the correlation plots
        """
        latex_file = open('paper/correlation_plots_si.txt', 'w')
        for i in range(len(self.variables_p)):
            latex_file.write(r'%')
            latex_file.write('\n')
            latex_file.write(r'\begin{figure}')
            latex_file.write('\n')
            latex_file.write(r'\includegraphics[width=')
            latex_file.write('{}'.format(0.45 / (7 - len(self.variables_p[i].split('_')))))
            latex_file.write(r'\textwidth]{post_' + self.variables_p[i] + '}')
            latex_file.write('\n')
            latex_file.write(r'\caption{\label{fig:post_' + self.variables_p[i] + '} ')
            latex_file.write(r'The posterior probabilities and correlation plots when ' + self.latex_p[i] + ' is varied in the model.}')
            latex_file.write('\n')
            latex_file.write(r'\end{figure}')
            latex_file.write('\n')
            latex_file.write(r'%')
            latex_file.write('\n')
        latex_file.close()

    def plot_refl_best(self):
        """
        Plot reflectometry when the best evidence.
        """
        file_h5 = h5py.File('output/model^{}.h5'.format(self.best_file), 'r')
        samples, weights = file_h5['samples'], np.exp(
            file_h5['logwt'] - file_h5['logz'][-1],
        )
        percentiles = np.zeros((samples.shape[1], 3))
        for i in range(samples.shape[1]):
            percentiles[i] = utils.quantile(
                samples[:, i],
                np.array([0.025, 0.5, 0.975]),
                weights,
            )

        variables = {}
        for i, var in enumerate(self.best_file.split('_')):
            variables[var] = i

        from model import refl

        axes = plt.subplots(figsize=(10, 15))[1]
        gobj = refl(
            percentiles[:, 1],
            CON_LIST,
            variables,
            False,
        )
        for i, obj in enumerate(gobj.objectives):
            axes.errorbar(
                obj.data.x,
                obj.data.y * obj.data.x ** 4 * 10 ** i,
                obj.data.y_err * obj.data.x ** 4 * 10 ** i,
                ls='',
                zorder=10,
            )
            axes.plot(
                obj.data.x,
                obj.model(obj.data.x) * obj.data.x ** 4 * 10 ** i,
                'k',
            )
        axes.set_yscale('log')
        axes.set_xlabel(r'$q$/Å$^{-1}$')
        axes.set_ylabel(r'$R(q)q^4$')
        plt.savefig('paper/refl.pdf')
        plt.close()

    def write_latex_percentiles(self):
        """
        Write latex snipets for the values from the best model.
        """
        file_h5 = h5py.File('output/model^{}.h5'.format(self.best_file), 'r')
        best_labels = self.latex_p[np.argsort(self.logz)[-1]]
        samples, weights = file_h5['samples'], np.exp(
            file_h5['logwt'] - file_h5['logz'][-1],
        )
        percentiles = np.zeros((samples.shape[1], 3))
        for i in range(samples.shape[1]):
            percentiles[i] = utils.quantile(
                samples[:, i],
                np.array([0.025, 0.5, 0.975]),
                weights,
            )
            file_latex = open(
                'paper/{}_range.txt'.format(
                    best_labels[1:-1].split('$/$')[i],
                ),
                'w',
            )
            file_latex.write(
                r'${:.1f}'.format(
                    percentiles[i, 1],
                ) + r'^{' + r'+{:.1f}'.format(
                    percentiles[i, 2] - percentiles[i, 1],
                ) + r'}_{' + r'-{:.1f}'.format(percentiles[i, 1] - percentiles[i, 0]) + r'}$',
            )
            file_latex.close()

plotter = Plotting()
plotter.generate_permutations()
plotter.make_plots_latex()
plotter.read_all_logz()
plotter.write_latex_logz()
plotter.write_latex_logz_table()
plotter.plot_logz()
plotter.plot_best_per_number()
plotter.plot_correlation_best()
plotter.plot_for_all()
plotter.plot_refl_best()
plotter.write_latex_percentiles()
for i in range(len(plotter.variables_p)):
    plotter.plot_for_any(plotter.variables_p[i], i)
