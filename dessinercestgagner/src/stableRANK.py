#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 20:38:10 2018

@author: Wojciech chacholski

Copyright Wojciech chacholski, 2018
This software is to be used only for activities related 
to SF2956 TDA course at KTH.
"""

import numbers
import numpy as np
inf=float("inf")

import matplotlib.pyplot as plt
plt.style.use('ggplot')

from ripser import ripser
###########################
## CLASSES
###########################

class pcf_max(object):
    """ The input:
        input_array: a 2xn array where:
            -the first row is the domain should be an increasing sequence of non-negative reals
            -the second row consists of values and any real value is allowed
            If the domain does not start with 0, the array will be augmented by 0,0"""
    def __init__(self,input_array):
        if np.asarray(input_array,dtype=float)[0][0]==0:
            self.content=np.asarray(input_array,dtype=float)
        else:
            self.content=np.insert(np.asarray(input_array,dtype=float),0,[0,0], axis=1)
    
    def pcf(self,precision=-2, base=10):
        """Converts a pcf_max into a pcf with a given precision and base"""
        if precision>=0:
            step=base**precision
        else:
            step=1/base**(-precision)
        last=np.int_(np.ceil(np.divide(self.content[0][-1],step)))
        D=np.arange(last+1)*step
        ind=np.searchsorted(self.content[0],D,side="right")-1
        return pcf(self.content[1][ind], precision,base)

    def simplify(self):
        """If the values are repeted for consequetive steps of the domain, then this step is removed"""
        C=self.content[1][:-1]-self.content[1][1:]
        k=np.insert(np.where(C!=0)[0]+1,0,0)
        outcome=self.content[:,k]
        return pcf_max(outcome)

    def evaluate(self,x):
        """evaluates self at x"""
        x_place=np.searchsorted(self.content[0], x, side='right')
        if  x<self.content[0][0]:
            return 0
        else:
            return self.content[1][x_place-1]    

    def dom_indices(self,b,e):
        '''gives the domain sequence between b and e. The sequence starts with b and ends with e'''
        if  b>e:
            raise ValueError("the indices can be calculated only for  [b e] in [0,\infty)")    
        else:
            indices=np.where(np.logical_and(self.content[0]>=b, self.content[0]<=e))[0]
            if len(indices)==0:
                dom=np.array([b,e])
            else:
                dom=self.content[0][indices]
                if self.content[0][indices[0]]>b:
                    dom=np.insert(dom, 0, b)
                if self.content[0][indices[-1]]<e:
                    dom=np.append(dom, e)
            return dom

    def restriction(self,b,e):
        '''restricts to the interval [b,e]'''
        indices=np.where(np.logical_and(self.content[0]>=b, self.content[0]<=e))[0] 
        if len(indices)==0:
            res=np.array([[b],[self.evaluate(b)]])
            #dom=np.array([b])
            #val=np.array([self.evaluate(b)])
        else:
            res=self.content[:,indices]
            #dom=self.content[0][indices]
            #val=self.content[1][indices]
            if res[0][0]>b:
                res=np.insert(res,0,[b,self.evaluate(b)],axis=1)
            #   res=np.insert(res,0,[b,self.evaluate(b)])
                #dom=np.insert(dom, 0, b)
                #val=np.insert(val, 0, self.evaluate(b))
            if res[0][-1]<e:
                res=np.append(res,[[e],[0]],axis=1)
                #dom=np.append(dom, e)
                #val=np.append(val, self.evaluate(e))
        return pcf_max(res)
    
    
    def extension(self,D):
        domain=np.unique(np.concatenate([self.content[0],D]))
        parameters=np.searchsorted(self.content[0],domain,side='right')-1
        val=self.content[1][parameters]
        return pcf_max(np.vstack( (domain, val) ))
          
    def integrate(self,b,e):
        ''' integrate  over the interval [b,e]'''
        integral=0
        l=0
        if b <= e:
            beg=b
            end=e
            coef=1
        else:
            beg=e  
            end=b 
            coef=-1
        temp = self.dom_indices(beg, end)
        integral = sum((temp[l+1]-temp[l])*self.evaluate(temp[l]) for l in range(len(temp) - 1))
        return coef*integral


    def plot(self,border=0,color="No",linewidth=1.0):   
        '''Plots with borders shifted to the left and right given by borders'''
        if border<0:
            raise ValueError("border has to be non-negative")
        x=self.content[0]
        y=self.content[1]
        if border!=0:
            x=np.append(x,x[-1]+border)
            y=np.append(y,y[-1])
        if color=="No":
            plt.plot(x, y, linewidth=linewidth, drawstyle='steps-post')
        else:
            plt.plot(x, y, color=color, linewidth=linewidth,drawstyle='steps-post')



    def __add__(self,other):
        '''adds two pcf_maxs or a number and a pcf_max'''
        if isinstance(other, numbers.Real):
            return pcf_max(np.vstack((self.content[0],self.content[1]+other)))        
        elif  type(other) == pcf_max:
            F1=self.extension(other.content[0])
            G1=other.extension(self.content[0])
            return  pcf_max(np.vstack((F1.content[0],F1.content[1]+G1.content[1])))
        else:
            raise ValueError("we can only add to pcf_max a pcf_max or a number")
        
    def __radd__(self,other):
        '''adds two pcf_maxs or number and pcf_max'''
        return self + other
        
    def __mul__(self, other):
        '''multily a pcf_max by a pcf_max or a number '''
        if isinstance(other, numbers.Real):
            return pcf_max(np.vstack((self.content[0],self.content[1]*other))) 
        elif  type(other) == pcf_max:
            F1=self.extension(other.content[0])
            G1=other.extension(self.content[0])
            return  pcf_max(np.vstack((F1.content[0],F1.content[1]*G1.content[1])))
        else:
            raise ValueError("we can only multiply pcf_max by a pcf_max or a number")            

       
    def __rmul__(self,other):
        '''multiply a pcf functions by  a number or a pcf function'''
        return self * other
        
    def __sub__(self, other):
        '''substracts for a pcf_max a pcf_max or a number'''
        if isinstance(other, numbers.Real):
            return pcf_max(np.vstack((self.content[0],self.content[1]-other)))        
        elif  type(other) == pcf_max:
            return  self+(-1)*other
        else:
            raise ValueError("we can substract from a pcf_max  a pcf_max  or a number")

    def invert(self):
        if np.any(self.content[1]==0):
            raise ValueError("we can not invert a pcf_max with a zero value")
        else:
            return pcf_max(np.vstack((self.content[0],1/self.content[1])))
            
    def __abs__(self):
        '''takes the absolute values of a pcf_max'''
        return pcf_max( np.vstack( (self.content[0], abs(self.content[1]) ) ))
    
    def apply(self,f):
        return pcf_max(np.vstack( (self.content[0], f(self.content[1]) ) ))
    
    def power(self, p):
        '''raises  pcf_max to a power'''
        return  pcf_max( np.vstack( (self.content[0], self.content[1] **p) ) ) 

    def lp_distance(self,other,p=1):
        ''' Gives an lp distance between two pcf_max functions'''
        if self.content[1][-1]!=other.content[1][-1]:
            return inf
        else:
            absdiff=abs(self - other)
            f=pcf_max.power(absdiff,p)
            b=absdiff.content[0][0]
            e=absdiff.content[0][-1]
            return (f.integrate(b,e))**(1/p)
        
        
    def interl(self,other):
        if np.all(self.content[1][1:] <= self.content[1][:-1])==False:
            raise ValueError("interleaving distance can be calculated only of pcf_maxs with  non increasing values")
        elif np.all(other.content[1][1:] <= other.content[1][:-1])==False:
            raise ValueError("interleaving distance can be calculated only of pcf_maxs with  non increasing  values")
        else:
            if self.content[1][-1]<other.content[1][-1]:
                return inf
            else:
                F1=self.extension(other.content[0])
                G1=other.extension(self.content[0])
                D=F1.content[0]
                F=F1.content[1]
                G=G1.content[1]
                FlG=np.where(F<G)[0]
                if len(FlG)==0:
                    return 0                
                else:
                    i=0
                    k=0
                    out=np.array([[FlG[0],FlG[0]]])
                    while i<len(FlG)-1:
                        if FlG[i+1]-FlG[i]==1:
                            out[k][1]+=1
                        else:
                            out=np.append(out,[[FlG[i+1],FlG[i+1]]],axis=0)
                            k+=1
                        i+=1

                    i=0
                    intervals=np.array([])
                
                    while i<len(out):
                        b=np.int(out[i][0])
                        e=np.int(out[i][1])+1
                        d=D[e]-D[b]
                        intervals=np.append(intervals,d)
                        i+=1
                        
                    return np.amax(intervals)

    def interleaving_distance(self,other):
        return max(pcf_max.interl(self,other),pcf_max.interl(other,self))

############################################################################
############################################################################
############################################################################
        
class pcf(object):
    """The input:
        input_vector could be either:
            - a vector of values, any real value is allowed
            - a 2xn array where:
                - the first row is the domain, should be an increasing sequence of non-negative reals
                - the second row consists of values and any real value is allowed
                If the domain does not start with 0, the array will be augmented by 0,0
        precision: it is an integer, the power to which the base will be raised,
            more negative more precise the calculations will be
        base: a positive integer, whose default is 10"""    

    def __init__(self, input_vector, precision=-2, base=10):
        if base<=0:
            raise ValueError("base should be a positive integer")
        if precision>=0:
            self.step=base**precision
        else:
            self.step=1/base**(-precision)
        self.precision=np.int_(precision)
        self.base=np.int_(base)     
        input=np.asarray(input_vector,dtype=float)        
        if input.ndim==2:
            if float(input_vector[0][0])!=0:
                input=np.insert(input,0,[0,0],axis=1)
            last=np.int_(np.ceil(np.divide(input[0][-1],self.step)))
            D=np.arange(last+1)*self.step
            ind=np.searchsorted(input[0],D,side="right")-1
            self.values=input[1][ind] 
            self.len=len(ind)        
        else:
            self.values=input
            self.len=len(input)


    def rescale(self, precision=-1):
        pr=np.int_(precision)
        if self.precision==pr:
            return self
        elif self.precision<pr:             
            coef=self.base**(pr-self.precision)
            values=np.mean(np.pad(self.values, (0, coef - self.len%coef), mode='constant', constant_values=self.values[-1]).reshape(-1, coef), axis=1)
            return pcf(values,precision=pr, base=self.base)
        else:
            coef=self.base**(self.precision-pr)
            values=self.values.repeat(coef)
            return pcf(values,precision=pr, base=self.base)        

    def pad(self,c):
        values=np.pad(self.values,(0,c),mode='constant', constant_values=self.values[-1])
        return pcf(values, precision=self.precision, base=self.base)    
        
    def reduct(self):
        """We can only reduct non-increasing PCFs"""
        if np.all(self.values[1:] <= self.values[:-1])==False:
            raise ValueError("We can only reduct non-increasing PCFs")
        k=np.where(self.values-self.values[-1] == 0)[0][0]
        outcome=self.values[0:k+1]
        return pcf(outcome,precision=self.precision, base=self.base)

    def plot(self,border=1, reduct="no",color="No",linewidth=1.0):
        if border<0:
            raise ValueError(" the border has to be  a non-negative number")
        if reduct=="yes":
            f=self.reduct()
        else:
            f=self
        x=np.arange(f.len)*f.step
        y=f.values
        if border!=0:
            x=np.append(x,x[-1]+border)
            y=np.append(y,y[-1])
        if color=="No":
            plt.plot(x, y, linewidth=linewidth, drawstyle='steps-post')
        else:
            plt.plot(x, y, color=color, linewidth=linewidth,drawstyle='steps-post')

    def integrate(self,b,e):
        step=self.step        
        if b==e:
            return 0
        else:
            if e>b:
                end=e
                beg=b
                coef=1
            else:
                end=b
                beg=e
                coef=-1                
            if end==inf:
                if self.values[-1]>0:
                    return coef*inf
                elif self.values[-1]<0:
                    return (-1)*coef*inf
                else:
                    end=(self.len-1)*step                    
            begdown=min(np.int_(np.floor_divide(beg,step)),self.len-1)
            enddown=min(np.int_(np.floor_divide(end,step)),self.len-1)            
            if begdown==enddown:
                return coef*(end-beg)*self.values[begdown]
            else:
                begup=np.int_(np.ceil(np.divide(beg,step)))
                I1=self.values[begdown]*(begup*step-beg)
                I3=self.values[enddown]*(end-enddown*step)
                I2=np.sum(self.values[begup:enddown])*step
                return coef*(I1+I2+ I3)
                


    def __add__(self,other):
        if isinstance(other, numbers.Real):
            return pcf(self.values+other,precision=self.precision, base=self.base)
        elif  type(other) == pcf: 
            if self.base!=other.base: 
                raise ValueError("we can add only pcf that have the same base")
            elif self.precision<other.precision:
                f=self
                g=other.rescale(self.precision)
                pr=self.precision
            elif self.precision>other.precision:
                g=other
                f=self.rescale(other.precision)
                pr=other.precision
            else:
                f=self
                g=other
                pr=self.precision
            if f.len<g.len:
                f=f.pad(g.len-f.len)
            elif f.len>g.len:
                g=g.pad(f.len-g.len)
            return  pcf(f.values+g.values,precision=pr,base=f.base)
        else:
            raise ValueError("we can add only  two pcf_reg or a pcf_reg and a number")
                
    def __radd__(self,other):
        '''adds two pcf functions or a number to a pcf function'''
        return self + other   

    def __mul__(self, other):
        if isinstance(other, numbers.Real):
            return pcf(self.values*other, precision=self.precision, base=self.base)
        elif  type(other) == pcf:
            if self.base!=other.base: 
                raise ValueError("we can muliply only pcf_reg that have the same base")
            elif self.precision<other.precision:
                f=self
                g=other.rescale(self.precision)
                pr=self.precision
            elif self.precision>other.precision:
                g=other
                f=self.rescale(other.precision)
                pr=other.precision
            else:
                f=self
                g=other
                pr=self.precision
            if f.len<g.len:
                f=f.pad(g.len-f.len)
            elif f.len>g.len:
                g=g.pad(f.len-g.len)
            return  pcf(f.values*g.values,precision=pr,base=f.base)
        else:
            raise ValueError("we can multiply to a pcf  a pcf_reg or a number")
       
    def __rmul__(self,other):
        '''multiply a pcf  by  a number or a pcf'''
        return self * other

    def __sub__(self,other):
        if isinstance(other, numbers.Real):
            return pcf(self.values-other,precision=self.precision, base=self.base)
        elif  type(other) == pcf: 
            if self.base!=other.base: 
                raise ValueError("we can subtract only pcf that have the same base")
            elif self.precision<other.precision:
                f=self
                g=other.rescale(self.precision)
                pr=self.precision
            elif self.precision>other.precision:
                g=other
                f=self.rescale(other.precision)
                pr=other.precision
            else:
                f=self
                g=other
                pr=self.precision
            if f.len<g.len:
                f=f.pad(g.len-f.len)
            elif f.len>g.len:
                g=g.pad(f.len-g.len)
            return  pcf(f.values-g.values,precision=pr,base=f.base)
        else:
            raise ValueError("we can subtract from a pcf only a pcf_reg or a number")

    def invert(self):
        if np.any(self.values==0):
            raise ValueError("we can invert a pcf with non-negative values")
        else:
            return pcf(1/self.values)
        
    def lp_distance(self,other,p=1):
        if self.base!=other.base: 
            raise ValueError("we can measure the lp distance of sig with the same base")
        elif abs(self.values[-1]-other.values[-1])>(1/self.base)**self.precision:   #self.values[-1]!=other.values[-1]:
            return inf
        elif self.precision<other.precision:
            f=self
            g=other.rescale(self.precision)
        elif self.precision>other.precision:
            g=other
            f=self.rescale(other.precision)
        else:
            f=self
            g=other
        if f.len<g.len:
            f=f.pad(g.len-f.len)
        elif f.len>g.len:
            g=g.pad(f.len-g.len) 
        return sum(f.step*np.absolute(f.values-g.values)**p)**(1/p)
        
    def interleaving_distance(self,other):
        if self.base!=other.base:
            raise ValueError("we can measure the interleaving distance of sig with the same base")
        if np.all(self.values[1:] <= self.values[:-1])==False:
            raise ValueError("interleaving distance can be calculated only of pcfs with  non increasing values")
        elif np.all(other.values[1:] <= other.values[:-1])==False:
            raise ValueError("interleaving distance can be calculated only of pcfs with  non increasing  values")
        else:
            if self.values[-1]!=other.values[-1]:
                return inf
            elif self.precision<other.precision:
                step=self.step
                f=self.values
                g=other.rescale(self.precision).values 
            elif self.precision>other.precision:
                step=other.step
                g=other.values
                f=self.rescale(other.precision).values
            else:
                step=self.step
                f=self.values
                g=other.values
            I=max(len(g)-np.searchsorted(g[::-1],f[::-1],side='right')[::-1]-np.arange(0,len(f)))
            J=max(len(f)-np.searchsorted(f[::-1],g[::-1],side='right')[::-1]-np.arange(0,len(g)))
            return max(I,J)*step


    def apply_function(self,f):
        return pcf(f(self.values),precision=self.precision, base=self.base)
    
    
############################################################################
############################################################################
############################################################################
        
class density(object):
    
    def __init__(self, input_array):
        input=np.asarray(input_array,dtype=float)
        if input.shape[0]!=2 or input.shape[1]==0:
            raise ValueError("the input for a density needs to be a 2xn array where n>0") 
        if np.any(input[1]<=0):
            raise ValueError("the values for a density have to be strictly positive")            
        if np.all(input[0][1:] > input[0][:-1])==False or input[0][0]<0:
            raise ValueError("the domain values for a density have to be strictly increasing and non-negative") 
        if input[0][0]>0:
            self.content=np.insert(input,0,[0,1], axis=1)
            self.dom=np.insert(input[0],0,0)
            self.values=np.insert(input[1],0,1)
        else:
            self.content=input
            self.dom=input[0]
            self.values=input[1]

    def evaluate(self,x):
        '''evaluates a density at x'''
        x_place=np.searchsorted(self.dom, x, side='right')
        if  x<self.dom[0]:
            return np.nan
        else:
            return self.values[x_place-1]       

    def restriction(self,b,e):
        '''restricts a density to the interval [b,e]'''
        indices=np.where(np.logical_and(self.dom>=b, self.dom<=e))[0] 
        if len(indices)==0:
            dom=np.array([b,e])
            val=np.array([self.evaluate(b),0])
        else:
            dom=self.dom[indices]
            val=self.values[indices]
            if self.dom[indices[0]]>b:
                dom=np.insert(dom, 0, b)
                val=np.insert(val, 0, self.evaluate(b))
            if self.dom[indices[-1]]<e:
                dom=np.append(dom, e)
                val=np.append(val, 0)
            else:
                val[-1]=0
        return np.vstack((dom,val))
            
    def integrate(self,b,e):
        ''' integrate a density function over the interval [b,e]'''
        if b <= e:
            beg=b
            end=e
            coef=1
        else:
            beg=e  
            end=b 
            coef=-1
        R=self.restriction(beg,end)
        return coef*sum(np.diff(R[0])*R[1][:-1])

    def reverse_integrate(self,a,t):
        if t<0:
            raise ValueError(" Reverse Integrate can only be done if the second parameter t is non negative")
        elif a>=self.dom[-1]:
            return a+t/(self.values[-1])
        else:
            f=self.restriction(a,self.dom[-1])
            CI=np.cumsum(np.diff(f[0])*f[1][:-1])            
            I=np.where(CI<=t)[0]
            if len(I)==0:
                return a+t/f[1][0]
            else:            
                ind=max(I)
                e=f[0][ind+1]
                return f[0][ind+1]+(t-CI[ind])/self.evaluate(e)

    def plot(self,border=1,color="No",linewidth=1.0):
        if border<0:
            raise ValueError(" the border has to be  a non-negative number")
        x=self.dom
        y=self.values

        if border!=0:
            x=np.append(x,x[-1]+border)
            y=np.append(y,y[-1])
        if color=="No":
            plt.plot(x, y, linewidth=linewidth, drawstyle='steps-post')
        else:
            plt.plot(x, y, color=color, linewidth=linewidth,drawstyle='steps-post')



    def dist_contour(self,a,t,truncate=(inf,inf)):
        trunc=truncate[0]
        trans_trunc=truncate[1]
        if a<0 or t<0 or trunc<=0 or trans_trunc<=0:
            raise ValueError("dis_contour is defined only on first quadrant and we can only trancate at positive numbers")
        else:
            RI=self.reverse_integrate(a,t)
            if RI<trunc and RI-a<trans_trunc:
                return RI
            else:
                return inf              

    def area_contour(self,a,t,truncate=(inf,inf)):
        trunc=truncate[0]
        trans_trunc=truncate[1]
        if a<0 or t<0 or trunc<=0 or trans_trunc<=0:
            raise ValueError("dis_contour is defined only on first quadrant and we can only trancate at positive numbers")
        else:
            y=self.reverse_integrate(0,a)
            C_at=a+self.integrate(y,y+t)
            if C_at<trunc and C_at - a<trans_trunc:
                return C_at
            else:
                return inf   

    def contour_plot(self,contour_type=("dist",inf,inf),scale=0.2):
        type=contour_type[0]
        truncate=contour_type[1:3]
        dom_end=max(self.dom[-1],self.reverse_integrate(0,1))
        int_end=self.integrate(0,dom_end)
        fig = plt.figure()
        self.plot()
        if type=="dist":
            a=np.arange(0,1.2*dom_end,scale)
            t=np.arange(0,1.2*int_end,scale)
            A, T = np.meshgrid(a,t)
            zs = np.array([self.dist_contour(a,t,truncate) for a,t in zip(np.ravel(A), np.ravel(T))])        
            Z = zs.reshape(A.shape)
            at=fig.add_subplot(111, projection='3d')
            at.plot_surface(A, T, Z)
            at.set_xlabel('a')
            at.set_ylabel('t')
            at.set_zlabel('dist contour')
        elif type=="area":
            a=np.arange(0,0.13*int_end,scale)
            t=np.arange(0,0.13*dom_end,scale)
            A, T = np.meshgrid(a,t)
            zs = np.array([self.area_contour(a,t,truncate) for a,t in zip(np.ravel(A), np.ravel(T))])        
            Z = zs.reshape(A.shape)
            at=fig.add_subplot(111, projection='3d')
            at.plot_surface(A, T, Z)
            at.set_xlabel('a')
            at.set_ylabel('t')
            at.set_zlabel('area contour')
        else:
             raise ValueError("you need to choose either dist or area for the contour type")                

    def bar_length(self,bar=np.array([1,2]), contour_type=("dist",inf,inf)):
        if len(bar)==0:
            return np.nan
        elif np.isnan(bar[1]) or bar[1]==inf:
            return inf
        else:
            birth=bar[0]
            death=bar[1]
            which=contour_type[0]
            trunc=contour_type[1]
            trans_trunc=contour_type[2]
            end=min([death,trunc,birth+trans_trunc])
            if birth>=end:
                return 0
            elif which=="dist":
                return self.integrate(birth,end)
            elif which=="area":
                return self.reverse_integrate(0,end)-self.reverse_integrate(0,birth)
            else:
                raise ValueError("you need to choose contour_type to be either dist or area ")   
                
                
                
###############################################################################
###############################################################################                
###############################################################################
                
class bc(object):
    
    def __init__(self,input):
        if isinstance(input,dict):
            self.input={}
            self.birth={}
            self.death={}
            for key in input:
                self.input[key]=np.asarray(input[key])
                if len(input[key])==0:
                    self.birth[key]=np.nan
                    self.death[key]=np.nan
                else:
                    self.birth[key]=np.asarray(input[key])[:,0]
                    self.death[key]=np.asarray(input[key])[:,1]
        else:
            self.input=np.asarray(input)
            if len(self.input)==0:
                self.birth=np.nan
                self.death=np.nan
            else:
                self.birth=self.input[:,0]
                self.death=self.input[:,1] 
            
            
            
    def length(self,contours={"H0": [ [[0],[1]], ("dist", inf, inf)] }):
        if isinstance(self.input,dict):
            outcome={}
            for key in self.input:
                if len(self.input[key])==0:
                    outcome[key]=np.nan
                else:
                    if key in contours.keys():
                        den=density(contours[key][0])
                        contour_type=contours[key][1]
                    else:
                        k=list(contours.keys())[-1]
                        den=density(contours[k][0])
                        contour_type=contours[k][1]
                    g= lambda x:  den.bar_length(x,contour_type)
                    outcome[key]= np.apply_along_axis(g,1,self.input[key])
            return outcome
        else:
            k=list(contours.keys())[0]
            den=density(contours[k][0])
            contour_type=contours[k][1]
            if len(self.input)==0:
                return np.nan
            else:        
                g= lambda x:  den.bar_length(x,contour_type)
                outcome=np.apply_along_axis(g,1,self.input)
            return outcome

    def stable_rank_max(self, contours={"H0": [  [[0],[1]], ("dist", inf, inf)] }):
        if isinstance(self.input,dict):
            length=self.length(contours)
            outcome={}
            for key in self.input:
                if len(self.input[key])==0:
                    outcome[key]=pcf_max(np.array([[0],[0]]))                    
                else:
                    sortlength=np.unique(length[key], return_counts=True)
                    dom=sortlength[0]
                    values=np.cumsum(sortlength[1][::-1])[::-1]
                    if dom[-1]==inf:
                        dom=np.insert(dom[:-1],0,0)
                    else:
                        dom=np.insert(dom,0,0)
                        values=np.append(values,0)   
                    outcome[key]=pcf_max(np.vstack((dom,values)))
            return outcome
        else:
            if len(self.input)==0:
                outcome=pcf_max(np.array([[0],[0]]))
            else:
                length=self.length(contours)
                sortlength=np.unique(length, return_counts=True)
                dom=sortlength[0]
                values=np.cumsum(sortlength[1][::-1])[::-1]
                if dom[-1]==inf:
                    dom=np.insert(dom[:-1],0,0)
                else:
                    dom=np.insert(dom,0,0) 
                    values=np.append(values,0)
                outcome=pcf_max(np.vstack((dom,values)))            
            return outcome




    
    def stable_rank(self, contours={"H0": [  [[0],[1]], ("dist", inf, inf)] }, precision=-2,base=10):
        return self.stable_rank_max(contours).pcf(precision,base)
        

        
    def plot(self, style="bar"):
        if isinstance(self.input,dict):
            if style=="bar":
                plt.yticks([])
                b=0
                for key in  self.input:
                    if np.isfinite(self.input[key]).all():
                        y=np.arange(b,b+len(self.input[key]))
                        plt.hlines(y,self.birth[key],self.death[key])
                        b1=y[-1]+np.int_(len(self.input[key])/20)
                        b=y[-1]+np.int_(len(self.input[key])/10)
                        lab=y[0]
                    else:
                        finite=self.input[key][np.isfinite(self.input[key]).all(axis=1)]
                        yfinite=np.arange(b,b+len(finite))
                        plt.hlines(yfinite,finite[:,0],finite[:,1])
                        a=np.amax(finite)
                        infinite=self.input[key][~np.isfinite(self.input[key]).all(axis=1)]
                        yinfinite=np.arange(b+len(finite),b+len(finite)+len(infinite))
                        plt.hlines(yinfinite,infinite[:0],np.full(len(infinite),a+1),color="blue")
                        b1=yinfinite[-1]+np.int_((len(finite)+len(infinite))/20)
                        b=yinfinite[-1]+np.int_((len(finite)+len(infinite))/10)
                        lab=yfinite[0]
                    plt.text(-1.4, lab, key, )
                    plt.axhline(b1+2, linestyle='dotted')
                    
            else:
                return print("for now only bar style")
        else:
            if style=="bar":
                y=np.arange(len(self.input))
                plt.hlines(y,self.birth,self.death)
            else:
                return print("for now only bar style")
                

#############################################################################
#############################################################################
#############################################################################
                            
class euc_object(object):
    
    def __init__(self,points=np.array([[0,0]]),  name=None, type=None):
        points=np.asarray(points,dtype=float)
        self.points=points
        self.number_points=len(points)
        self.dim=points.shape[1]
        self.name=name
        self.type=type
        
    def plot(self, color="No"):      
        dim=self.dim
        if dim==2:     
            X=self.points[:,0]
            Y=self.points[:,1]
            if color=="No":
                return plt.scatter(X,Y)
            else:
                return plt.scatter(X,Y,c=color)
        elif dim==3:
            fig = plt.figure(self.name)
            ax = fig.add_subplot(111, projection='3d')
            X=self.points[:,0]
            Y=self.points[:,1]
            Z=self.points[:,2]
            if color=="No":    
                return ax.scatter(X,Y,Z)
            else:
                return ax.scatter(X,Y,Z,c=color)
        else:
            raise ValueError("We can only plot 2 and 3 d objects")
     
    def distance_matrix(self, metric="euclidean"):
        return spatial.distance.pdist(self.points,metric)
    


    def barcode(self, maxdim=1, thresh=inf, coeff=2, distance_matrix=False, do_cocycles=False, metric='euclidean'):
        dgms=ripser(self.points, maxdim, thresh, coeff, distance_matrix, do_cocycles, metric)["dgms"]
        i=0
        dic={}
        while i<=maxdim:
            dic["H"+str(i)]=dgms[i]
            i=i+1
        return bc(dic)    

    def union(self,other):
        if self.dim==other.dim:
            spoints=self.points
            opoints=other.points

        elif self.dim<other.dim:
            d=other.dim-self.dim
            z=np.zeros((self.number_points, d), dtype=self.points.dtype)            
            spoints=np.concatenate((self.points,z), axis=1)
            opoints=other.points
        else:
            d=self.dim-other.dim
            z=np.zeros((other.number_points, d), dtype=self.points.dtype)
            opoints=np.concatenate((other.points,z), axis=1)
            spoints=self.points
        
        allpoints=np.concatenate((spoints,opoints))
        points=np.unique(allpoints,axis=0)    
        return euc_object(points)
                
    
    def translate(self,vec):
        if len(vec)>= len(self.points):
            return euc_object(self.points+vec[:len(self.points)])
        else:
            vector=np.pad(vec,(0,self.dim-len(vec)),"constant",constant_values=0)
            return euc_object(self.points+vector)

    def sub_sample(self,number_points):
        if number_points>=self.number_points:
            return self
        else:
            A=np.random.choice(np.arange(self.number_points), number_points, replace=False)
            points=self.points[A]
            return euc_object(points,  self.name, self.type)
            


    


    