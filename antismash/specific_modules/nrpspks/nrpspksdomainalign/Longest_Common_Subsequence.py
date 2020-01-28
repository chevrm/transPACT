#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Longest Common Subsequence, Shift algorithm.
# Copyright (C) 2012  Gonzalo Exequiel Pedone
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with This program. If not, see <http://www.gnu.org/licenses/>.
#
# Email   : hipersayan DOT x AT gmail DOT com
# Web-Site: http://hipersayanx.blogspot.com/
#
# The shift algorithm consist on finding the Longest Common Subsequence between
# a sequence 'a' and 'b' by comparing the last element of 'a' with the first
# element of 'b', and finding the longest sequence between the common elements,
# then shifting 'a' to the right, and comparing the next two elements, and so
# on, until reach the right wall, then shifting 'b' to the left until reach the
# left wall.
# This is a non-recursive algorithm, and the maximum number of comparisons are:
# (len(a) + len(b) - 1)**2.
#
#                 ===>
#             __ __ __ __
#          a |__|__|__|__|__ __ ) right wall
# left wall (         |__|__|__| b
#                        <===

def lcs(a=[], b=[]):
    if a == [] or b == []:
        return []

    l = len(a) + len(b) - 1

    # Fill non-comparable elements with null spaces.
    sa = a + (len(b) - 1) * ['']
    sb = (len(a) - 1) * [''] + b

    longest = []

    for k in range(l):
        cur = []

        for c in range(l):
            if sa[c] != '' and sb[c] != '' and sa[c] == sb[c]:
                cur.append(sa[c])
            else:
                if len(cur) > len(longest):
                    longest = cur

                cur = []

        if len(cur) > len(longest):
            longest = cur

        if sa[len(sa) - 1] == '':
            # Shift 'a' to the right.
            sa = [''] + sa[: len(sa) - 1]
        else:
            # Shift 'b' to the left.
            sb = sb[1:] + ['']

    return longest

'''
a = 'Andrea always wears black cat to all places she go.'
b = 'She likes her black cat because it is very mischievous and playful.'

print(a)
print(b)
print()
print('What carry Andrea always?')
print()
print(''.join(lcs(list(a), list(b))))
'''




# This sequences alignment algorithm consist in finding the LCS of two sequences
# and then resolve the alignment of the sequences before and after the LCS
# recursively.

# Auxiliary function for finding the first index in which appears a sub-sequence.
def findSubList(l=[], sub=[]):
    if len(sub) > len(l):
        return -1

    for i in range(len(l) - len(sub) + 1):
        j = 0
        eq = True

        for s in sub:
            if l[i + j] != s:
                eq = False

                break

            j += 1

        if eq:
            return i

    return -1

def alignSequences(sequence1=[], sequence2=[]):
    # lcs is the Longest Common Subsequence function.
    cs = lcs(sequence1, sequence2)

    if cs == []:
        return sequence1 + [''] * len(sequence2) , \
               [''] * len(sequence1) + sequence2
    else:
        # Remainings non-aligned sequences in the left side.
        left1 = sequence1[: findSubList(sequence1, cs)]
        left2 = sequence2[: findSubList(sequence2, cs)]

        # Remainings non-aligned sequences in the right side.
        right1 = sequence1[findSubList(sequence1, cs) + len(cs):]
        right2 = sequence2[findSubList(sequence2, cs) + len(cs):]

        # Align the sequences in the left and right sides:
        leftAlign = alignSequences(left1, left2)
        rightAlign = alignSequences(right1, right2)

        return leftAlign[0] + cs + rightAlign[0], \
               leftAlign[1] + cs + rightAlign[1]

'''
X = ['IV', 'III_b', 'III_e', 'I_a','XI_b', 'IX', 'II_d', 'II_c3', 'VI_c',	'XII_a', 'III_a']
Y = ['IV_e', 'III_b', 'III_e', 'I_a', 'XI_b', 'IX', 'II_d',	'II_c3', 'XII_a', 'III_a','noInfo']
a, b = alignSequences(list('abcdfghjqz'), list('abcdefgijkrxyz'))
a, b = alignSequences(Y, X)
print(':'.join(['-' if i == '' else i for i in a]))
print(':'.join(['-' if j == '' else j for j in b]))
'''






if __name__ == "__main__":
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    parser = argparse.ArgumentParser()
    parser.add_argument('--X')
    parser.add_argument('--Y')

    args = parser.parse_args()
    if any(t == [] for t in vars(args).values()) == False:
        X = args.X
        Y = args.Y
        a, b = alignSequences(Y, X)
        print(':'.join(['-' if i == '' else i for i in a]))
        print(':'.join(['-' if j == '' else j for j in b]))
    else:
        print """Please provide needed profiles"""
        sys.exit()
