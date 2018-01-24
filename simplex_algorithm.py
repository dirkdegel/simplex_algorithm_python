# status: solves feasible and bounded (linear) problems in canonical form: max{c^T x | Ax <= b, x >= 0} with b >= 0
# tested: Python 3.6.1

# to do's
# major:
# 1) unbounded check -> how to implement? in each case or only if column is selected as pivot column?
# 2) primal and dual degeneration check
# 3) switch to first improvement strategy in case of stalling
# x) two phase simplex
# minor:
# theta calculation in case of lhs = 0
# display only 2 digits in display()
# check if minimization works correct

from __future__ import division
from numpy import *


class tableau:

    def __init__(self):
        self.obj = []
        self.rows = []
        self.cons = []

    def add_objective(self, objective_function, sense):
        self.obj.append(sense)
        self.obj += objective_function

    def add_constraint(self, constraint, rhs):
        self.rows.append([0] + constraint)
        self.cons.append(rhs)

    def construct(self):
        # signs of objective function depending on minimization (-1) or maximization (+1)
        if self.obj[0] == 1:
            for j in range(1,len(self.obj)):
                self.obj[j] *= -1
        for i in range(len(self.rows)):
            self.obj += [0]
            unit_vector = [0 for k in range(len(self.rows))]
            unit_vector[i] = 1
            self.rows[i] += unit_vector + [self.cons[i]]
            self.rows[i] = array(self.rows[i], dtype=float)
        self.obj += [0]

    def display(self):
        print('\n', matrix([self.obj] + self.rows))

    # determine pivot column (steepest edge rule)
    def _calc_pivot_column(self):
        delta_z = 0
        pivot_column = 0
        for j in range(1, len(self.obj) - 1):   # [0,...,n-1]
            if self.obj[j] < delta_z:
                delta_z = self.obj[j]
                pivot_column = j
        #if pivot_column == 0: return -1     # tableau is optimal # is not needed do to _optimality_check()
        is_bounded = 0
        for i in range(len(self.rows)):
            if self.rows[i][pivot_column] > 0:
                is_bounded = 1
            #else:
            #    continue
        if is_bounded == 1:
            return pivot_column
        else:
            return -1

    # determine pivot row
    # bounded check could be included here?!
    def _calc_pivot_row(self, pivot_column):
        rhs = [self.rows[i][-1] for i in range(len(self.rows))]     # -1 is the last element in the row, i.e. the RHS
        lhs = [self.rows[i][pivot_column] for i in range(len(self.rows))]
        theta = []
        for i in range(len(rhs)):   # [0,...,m-1]
            if lhs[i] > 0:
                theta.append(rhs[i] / lhs[i])
            else:
                theta.append(99999999 * abs(max(rhs)))  # find better formulation
                #theta.append(max(rhs) / min(abs(lhs)))  # check it !!!
        return argmin(theta)    # pivot row

    # calculate pivot operations
    def _pivot(self, pivot_row, pivot_column):
        pivot_element = self.rows[pivot_row][pivot_column]
        self.rows[pivot_row] /= pivot_element
        for r in range(len(self.rows)):
            if r == pivot_row:
                continue
            else:
                self.rows[r] = self.rows[r] - self.rows[r][pivot_column] * self.rows[pivot_row]
        self.obj = self.obj - self.obj[pivot_column] * self.rows[pivot_row]

    # optimality check
    def _optimality_check(self):
        if min(self.obj[1:-1]) >= 0:
            return 1   # solution is optimal
        else:
            return 0   # perform next simplex iteration

    def _dual_degeneration_check(self):
        if count_nonzero(self.obj[1:-1]) < len(self.obj) - 2 - len(self.rows):
            return 1    # is dual degenerated
        #else:
            #return 0    # is not dual degenerated

    def _primal_degeneration_check(self):
        if count_nonzero([self.rows[i][-1] for i in range(len(self.rows))]) < len(self.rows):
            return 1    # is primal degenerated

    # simplex algorithm
    def solve(self):
        self.display()  # initial tableau
        while not self._optimality_check():
            c = self._calc_pivot_column()
            if c == -1:
                print('problem is unbounded')
                break
            r = self._calc_pivot_row(c)
            self._pivot(r, c)
            #print('\npivot column: %s\npivot row: %s' % (c + 1, r + 2)) # figure out what is more intuitive
            print('\npivot column: %s\npivot row: %s' % (c, r))
            self.display()
            if self._primal_degeneration_check() == 1:
                print('problem is primal degenerated!')
        if self._dual_degeneration_check() == 1:
            print('problem is dual degenerated!')

# --- end of the simplex algorithm ---




# problem input

# standard case
t = tableau()
t.add_objective([2, 3, 2], 1)
t.add_constraint([2, 1, 1], 4)
t.add_constraint([1, 2, 1], 7)
t.add_constraint([0, 0, 1], 5)
t.construct()
#t.solve()

# standard case (all edges get visited, Klee–Minty cube)
t2 = tableau()
t2.add_objective([2, 1], 1)
t2.add_constraint([1, 0], 5)
t2.add_constraint([4, 1], 25)
t2.construct()
#t2.solve()

# unbounded, solvable

# unbounded from Werners, Brigitte - Grundlagen des Operations Research
t3 = tableau()
t3.add_objective([2, 4], 1)
t3.add_constraint([-2, 3], 12)
t3.add_constraint([-1, 3], 18)
t3.construct()
#t3.solve()

# unfeasible

# dual degeneration from Werners, Brigitte - Grundlagen des Operations Research
t4 = tableau()
t4.add_objective([3, 1.5], 1)
t4.add_constraint([2, 1], 22)
t4.add_constraint([1, 2], 23)
t4.add_constraint([4, 1], 40)
t4.construct()
#t4.solve()

# primal degeneration from Werners, Brigitte - Grundlagen des Operations Research
t5 = tableau()
t5.add_objective([3, 2], 1)
t5.add_constraint([2, 1], 22)
t5.add_constraint([1, 2], 23)
t5.add_constraint([4, 1], 40)
t5.add_constraint([2, 0.75], 21)
t5.construct()
t5.solve()

'''

if __name__ == '__main__':

    t = Tableau()
    t.add_objective([2, 3, 2], 1)
    t.add_constraint([2, 1, 1], 4)
    t.add_constraint([1, 2, 1], 7)
    t.add_constraint([0, 0, 1], 5)
    t.construct()
    t.solve()

'''