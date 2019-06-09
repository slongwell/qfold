#!/usr/bin/env python
# coding: utf-8

# This notebook is converted into a module with the following command:   
# `jupyter nbconvert --to script hamiltonians.ipynb`

# # General

# In[1]:


import numpy as np
from pyquil.paulis import sI, sX, sY, sZ

global N
global D
global n
global q

# setter for global variables
def init_ham(N_,D_):
    global N
    global D
    global n
    global q
    N = N_
    D = D_
    if D_ == 2:
        n = 4*N_-10
        q = [sZ(x) for x in range(n)]
        
# inclusive range (for sums, products)
def inc(start,stop):
    return range(int(start),int(stop+1))


# $\begin{array}{rl}
# N & \text{# residues in protein}\\
# n & \text{# qubits}\\
# D & \text{# dimensions}\\
# \end{array}$
# 
# 0-indexed, start inclusive, stop exclusive:
# 
# $\begin{array}{rl}
# i,i' \in [0,N)  & \text{protein residue indices}\\
# j,j' \in [0,n) & \text{qubit indices}\\
# k,k' \in [0,2D) & \text{direction indices}\\
# t \in [0,N-1)    & \text{turn index}\\
# \end{array}$
# 
# $\begin{array}{rl}
# P_{i,i'} & \text{contact energy between adjacent residues $i$ and $i'$}
# \end{array}$

# In[22]:


def AND(x,y):
    return x*y


# $\texttt{AND}(x,y) = xy$   

# In[25]:


def XNOR(x,y):
    global q
    return 1 - x - y + 2*x*y


# $\texttt{XNOR}(x,y) = 1 -x - y + 2xy$

# In[38]:


def XOR(x,y):
    return x + y - 2*x*y


# $\texttt{XOR}(x,y) = x + y - 2xy$

# In[41]:


def HA(x,y):
    return AND(x,y), XOR(x,y)


# $\texttt{HA}(x,y) = (xy,\text{ } x+y-2xy)$

# In[217]:


def op(k):
    global D
    return (k+D)%(2*D)

def d(k,t):
    global q
    if (t==1) & (k in [0,1]):
        return q[k]
    elif (t==1) & (k not in [0,1]):
        return 0
    else:
        return q[4+k-6]


# bitstrings: left-indexed
# 
# $k_o = (k+D) \text{ mod } (2D)$
# 
# 2D:
# 
# $d_k^t = 
# \cases{
# q_k        & $t=1, k \in \{0,1\}$ \cr
# 0          & $t=1, k \not\in \{0,1\}$ \cr
# q_{4t+k-6} & else
# }$

# In[ ]:


# return sum string from list of bits
# this circuit does not include the quasilinear optimization
def S(lst, verbose=False): 
    len_out = int(np.ceil(np.log2(len(lst)+1)))
    for i in range(len(lst)-1):
        start = len(lst) - i%2
        stop = start-i-1
        for j in range(start,stop,-2):
            if verbose:
                print(i*' ', j,j-1)
            lst[j-1],lst[j-2] = HA(lst[j-1],lst[j-2])
            
    for i in reversed(range(len(lst)-2)):
        start = len(lst) - i%2
        stop = start-i-1
        for j in range(start,stop,-2):
            if verbose:
                print(i*' ', j,j-1)
            lst[j-1],lst[j-2] = HA(lst[j-1],lst[j-2])
    return lst[:len_out]


# ![adder](images/adder.png)

# In[2]:


def s(k,r,i,ii):
    lst = [d(k,t) for t in range(4)]
    return S(lst)[r-1]


# $s_k^r(i,i') = r^{th} \text{digit of} \sum_{p=i}^{i'-1} d_k^p$
# 
# $\begin{array}{rl}
# r & \text{sumstring index}\\
# p & \text{'rightmost' 0 in the odd, lesser sum string}\\
# \end{array}$
# - sum strings (*Babbush2013: 2.3.1*): 
#     - ancilla qubits?
#     - $1 \le i \lt i' \le N$
#     - len $log_2(i' - i)$

# -----------
# -----------
# ## $H_{overlap}$

# In[155]:


def h_overlap(i,ii):
    global D
    return np.prod([XNOR(s(k,r,i,ii),s(op(k),r,i,ii)) 
                   for k in inc(1,D) for r in inc(1,np.ceil(np.log2(ii-i)))])


# $h_{overlap}(i,i') = \prod_{k=1}^{D} \prod_{r=1}^{\lceil log_2(i'-i) \rceil} \texttt{XNOR}(s_k^r(i,i'),s_{k_o}^r(i,i'))$

# In[137]:


def H_overlap(lambda_overlap):
    global N
    return lambda_overlap * sum(h_overlap(i,i+2*ii) for i in inc(0,N-3) for ii in inc(1,np.floor(0.5*(N-i-1))))


# $\boxed{
# H_{overlap} = \lambda_{overlap} \sum_{i=0}^{N-3} \sum_{i'=1}^{\lfloor \frac{1}{2}(N-i-1) \rfloor} h_{overlap}(i,i+2i')
# }$
# - To set $\lambda_{overlap}$, see [*Babbush2013: 7.2.4*](https://arxiv.org/pdf/1211.3422.pdf)

# ## $H_{pair}$

# In[ ]:


def t120(k,i,ii,p):
    return XOR(s(k,p-1,i,ii),s(k,p,i,ii))


# $t_{120}(k,i,i',p) = \texttt{XOR}(s_{k}^{p-1}(i,i'),s_{k}^p(i,i'))$

# In[ ]:


def t121(k,i,ii,p):
    return np.prod([XNOR(s(k,r,i,ii),s(k,r+1,i,ii))
                    for r in inc(1,p-2)])


# $t_{121}(k,i,i',p) = \prod_{r=1}^{p-2} \texttt{XNOR}(s_{k}^r(i,i'),s_{k}^{r+1}(i,i'))$

# In[ ]:


def t122(k,i,ii,p):
    return np.prod([XOR(s(k,r,i,ii),s(op(k),r,i,ii))
                    for r in inc(1,p)])


# $t_{122}(k,i,i',p) = \prod_{r=1}^p \texttt{XOR}(s_{k}^r(i,i'),s_{k_{o}}^r(i,i'))$

# In[ ]:


def t123(k,i,ii,p):
    return np.prod([XNOR(s(k,r,i,ii),s(op(k),r,i,ii))
                    for r in inc(p+1,np.ceil(np.log2(ii-i)))])


# $t_{123}(k,i,i',p) = \prod_{r=p+1}^{\lceil log_2(i'-i) \rceil} \texttt{XNOR}(s_{k}^r(i,i'),s_{k_{o}}^r(i,i'))$

# -----

# In[ ]:


def t10(k,i,ii):
    return XOR(s(k,1,i,ii),s(op(k),1,i,ii))


# $t_{10}(k,i,i') = \texttt{XOR}(s_{k}^1(i,i'),s_{k_o}^1(i,i'))$

# In[ ]:


def t11(k,i,ii):
    return np.prod([XNOR(s(k,r,i,ii),s(op(k),r,i,ii))
                    for r in inc(2,np.ceil(np.log2(ii-i)))])


# $t_{11}(k,i,i') = \prod_{r=2}^{\lceil log_2(i'-i) \rceil}\texttt{XNOR}(s_{k}^r(i,i'),s_{k_o}^r(i,i'))$

# In[ ]:


def t12(k,i,ii):
    return sum(t120(k,i,ii,p) * t121(k,i,ii,p) * t122(k,i,ii,p) * t123(k,i,ii,p) for p in inc(2,np.ceil(np.log2(ii-i)))) 


# $t_{12}(k,i,i') = \sum_{p=2}^{\lceil log_2(i'-i) \rceil} t_{120}(k,i,i') * t_{121}(k,i,i')* t_{122}(k,i,i')* t_{123}(k,i,i')$

# ----

# In[ ]:


def t0(k,i,ii):
    kk_lst = [kk for kk in range(2*D) if kk != k]
    return np.prod([XNOR(s(kk,r,i,ii),s(op(kk),r,i,ii)) 
                   for kk in kk_lst for r in inc(1,np.ceil(np.log2(ii-i)))])


# $t_0(k,i,i') = \prod_{k'\neq k} \prod_{r=1}^{\lceil log_2(i'-i) \rceil} \texttt{XNOR}(s_{k'}^r(i,i'),s_{k_o'}^r(i,i'))$

# In[ ]:


def t1(k,i,ii):
    return t10(k,i,ii) * t11(k,i,ii) + t12(k,i,ii)


# $t_1(k,i,i') = t_{10}(k,i,i') * t_{11}(k,i,i') + t_{12}(k,i,i')$

# ----

# In[ ]:


def a(k,i,ii):
    return t0(k,i,ii)*t1(k,i,ii)


# $a_k(i,i') = t_0(k,i,i') * t_1(k,i,i')$
# 
# ![](images/adjacency.png)

# ----

# In[222]:


def H_pair(P):
    global N
    global D
    return sum(P[i -1,1+i+2*ii -1] * sum(a(k,i,1+i+2*ii) for k in inc(1,D)) for i in inc(1,N-3) for ii in inc(1,np.floor(0.5*(N-i-1))))


# $\boxed{H_{pair} = \sum_{i=1}^{N-3} \sum_{i'=1}^{\lfloor \frac{1}{2}(N-i-1)\rfloor} P_{i,1+i+2i'} \sum_{k=1}^{D} a_k(i,1+i+2i')}$
# - [*Fingerhuth 2018*](https://arxiv.org/pdf/1810.13411.pdf) is missing floor, has lower case $N$ in first sum
# - [*Babbush 2013*, eqn (40)](https://arxiv.org/pdf/1211.3422.pdf)

# # $H_{cost}$

# In[3]:


# YY
def H_cost(lambda_overlap,P):
    return H_overlap(lambda_overlap) + H_pair(P)


# $\boxed{H_{cost} = H_{overlap} + H_{pair}}$

# -----------
# -----------
# # $H_{mixer}$ Variants:

# In[1]:


def X_mixer():
    global n
    return sum(sX(i) for i in range(0,n))


# $\boxed{X_{mixer} = \sum_{i} X_i}$

# ----

# In[10]:


def XY(i,ii):
    return 0.5*(sX(i)*sX(ii) + sY(i)*sY(ii))


# $\texttt{XY}_{i,i'} = \frac{1}{2}(X_i X_{i'} + Y_i Y_{i'})$
# - has similar properties to a $\texttt{SWAP}$ gate (actually redefined as $\texttt{SWAP}$ in *Fingerhuth 2018*)

# In[ ]:


def M(t):
    global n
    return sum(XY(i,ii) for i in inc(0,n-2) for ii in inc(i+1, n-1))


# $M_t = \sum_{i=0}^{n-2} \sum_{i'=i+1}^{n-1} \texttt{XY}_{i,i'}$

# In[1]:


def XY_simple():
    global N
    return sum(M(t) for t in inc(1,N-2))


# $\boxed{XY_{simple} = \sum_{t=1}^{N-2} M_t}$
