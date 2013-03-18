# https://github.com/mozilla/datazilla-metrics/blob/master/dzmetrics/fdr.py
# SHA 3168305e766e109f5ad2d6690aad7b6261b5931f
#
# Author: Joseph Kelly, Mozilla metrics.
#
# License (https://github.com/mozilla/datazilla-metrics/blob/master/LICENSE.txt)
#
# Copyright (c) 2012 Mozilla
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of the author nor the names of other
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

def fdr(p_values, q=0.05):
    """
    Implements the Benjamini-Hochberg method of false discovery rate control.

    See e.g. http://en.wikipedia.org/wiki/False_discovery_rate,
    http://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg_procedure
    and http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va

    Given a list of p-values (floats) for independent comparisons, and a q
    value (the upper bound on the false discovery rate; the expected proportion
    of false rejections of the null hypothesis), returns a dictionary with two
    keys: "status" is a list of boolean values the same length as the given
    list of p-values, where a ``True`` value represents rejection of the null
    hypothesis for that p-value and ``False`` represents acceptance of the null
    hypothesis, and "count" is the number of p-values for which the null
    hypothesis was rejected (these will always be the lowest "count" p-values).

    Same as p.adjust(pvs, method = 'BH') but returns an array of bools: True if hypothesis is rejected
    """

    N = len(p_values)
    index = range(0, N)
    #pindex = zip(p_values, index)
    #sortedp = sorted(pindex)
    sortedp = sorted(zip(p_values, index))
    #print "DEBUG: sortedp = %s" % sortedp

    # find cutoff for rejection
    cutoff = [(i+1)*q/float(N) for i in index]
    indicator = 0
    for i in index:
        if sortedp[i][0] < cutoff[i]:
            indicator = i + 1
    #print "DEBUG: indicator = %s" % indicator

    # reject/fail to reject
    status = [True]*indicator + [False]*(N-indicator)
    output = range(0, N)
    for i in index:
        output[sortedp[i][1]] = status[i]

    #return {"status": output, "count": indicator}
    return output

