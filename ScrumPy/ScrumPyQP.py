"""
QP implementation by Thomas Pfau based on ScrumPyLP by Mark Poolman
"""

import sys, types, exceptions


SysStderr = sys.stderr
NullF =  open("/dev/null", "w")


#Some definitions of static variables
GREATER_THAN = "G"
LESS_THAN = "L"
EQUAL_TO = "E"
RANGED = "R"
import cplex
import math
from Util import Set, DynMatrix
from Data import DataSets
from Structural import StoMat
from ScrumPy import Model as ModelClass

########################### Mappings for human readability
Status2Str = {           # map symbolic status constants to human strings
    1    : "Optimal",
    23   : "Feasible",
    3 : "Infeasible",
    5 : "No feasible",
    2  : "Unbounded!",
    4  : "Unbounded!",
    101 : "Optimal",
    102 : "Optimal", #within tolerance
    103 : "Infeasible",
    }

ObjectiveMap = {
    "Max" : cplex.Cplex.objective.sense.maximize,
    "Min" : cplex.Cplex.objective.sense.minimize,
    cplex.Cplex.objective.sense.minimize : "Min" ,
    cplex.Cplex.objective.sense.minimize : "Max"
    }
###############################

VerySmall = 1e-7
def IsZero(number, tol=VerySmall):
    return abs(number) < tol

###############################

def SplitRev(sm,Conv=float): # copy of sm with reversible reacs split into forward and back irreversible pairs

    rv = sm.Copy(Conv)
    rv.Reversed = reversed = []
    rv.Backs = backs = {}
    
    for r in sm.cnames:
#        print "Setting up Reaction " + r
        if sm.RevProps[r] == StoMat.t_Rever:
            sto = sm.InvolvedWith(r)
            rback = r+"_back"
            try:
                rv.NewReaction(rback,sto)
            except:
                print r
                raise "Error creating new Reaction"
            rv.MulCol(rback,k=-1)
            backs[r] = rback

        elif sm.RevProps[r] == StoMat.t_BackIrrev:
           sm.MulCol(reac,  k=-1)
           reversed.append(reac)

    return rv


class qp(cplex.Cplex):

    def __init__(self, model,infinity = cplex.infinity):
        global VerySmall
        cplex.Cplex.__init__(self)
        self.infinity = infinity
        self.updateActive = True;
        self.set_problem_name("Quadratic Problem")
        self.set_results_stream(sys.__stdout__)
        if not isinstance(model,ModelClass):
            self.variables_names = []
            self.linear_constraints_names = []
            return
        self.sm = SplitRev(model.sm)
        #cplex.Cplex.__init__(self)
        #Set all reactions to be greater than zero        
        self.variables.add(names=self.sm.cnames[:],lb=len(self.sm.cnames)*[0],ub = len(self.sm.cnames)*[self.infinity])
        self.variables_names = self.sm.cnames[:] # print self.variables.get_names()
        
        #setup the rows (i.e. constraints)
        self.linear_constraints.add(names=self.sm.rnames[:])
        self.linear_constraints_names = self.sm.rnames[:]
        self.quadratic_constraints_names = []
        for row in self.sm.rnames:
  #          print "Setting up Row: " + row
            self.SetRowVals(row, map(float, self.sm[row]))
            self.SetRowBounds(row, 0,0)
        self.rnames = {}
        self.ridxs = {}
        self.cnames = {}
        self.cidxs = {}
        self.updateRowDics()
        self.updateColDics()
        VerySmall = self.parameters.feasopt.tolerance.get()

        
    def __del__(self):
        cplex.Cplex.__del__(self)

    def __getCol(self, c):
        if c in self.variables_names:
            return self.variables_names.index(c)
        else:
            raise exceptions.KeyError, str(c)
        
    def __getRow(self, r):
        if r in self.linear_constraints_names:
            return self.linear_constraints_names.index(c)
        else:
            raise exceptions.KeyError, str(c)
      
    #compatability functions, these are not really necessary
    def updateColDics(self):
        if not self.updateActive:
            return
        self.cnames = dict(zip(range(len(self.variables_names)),self.variables_names))
        self.cidxs =  dict(zip(self.variables_names,range(len(self.variables_names))))

    def updateRowDics(self):
        if not self.updateActive:
            return
        """Be aware, that this function is only for linear constraints!!"""
        self.rnames = dict(zip(range(len(self.linear_constraints_names)),self.linear_constraints_names))
        self.ridxs = dict(zip(self.linear_constraints_names,range(len(self.linear_constraints_names))))
    
    ##
    ##
    #   Objective Modifications:
    ##
    ##

    ##
    #   General
    ##

    def AddToObjective(self,reacs,quadratic=False):
        if quadratic:
            self.SetQuadObjective(reacs)
        else:
            self.SetLinObjective(reacs)
    
    def SetObjective(self,reacs,quadratic=False):
        self.ClearObjective()
        if quadratic:
            self.SetQuadObjective(reacs)
        else:
            self.SetLinObjective(reacs)

    def createDiffObjectiveForReacs(self,reacs):
        """ reacs has to be a dictionary of reaction fluxes (w). all other fluxes are assumed to be 0"""
        self.ClearObjective()
        linpart = reacs.copy()
        #the direction is minimising!
        self.SetObjDirec("Min")
        for reac in reacs:
            back = self._BackReac(reac)
            if back != None:
                linpart[back] = linpart[reac]
            linpart[reac] = -linpart[reac]
            if not reac in self.variables_names:
                print reac + " not in the problem"
                continue
        #This means we get -w*x + x^2 as objective. which is, if minimized the same as minimizing the eukleadian norm.
        self.SetLinObjCoefsFromDic(linpart)
        self.SetQuadObjCoefsFromLists(reacs.keys(),[0.5]*len(reacs.keys()))

    def CreateIndicatorConstraintsAndObjective(self,reacs):
        self.ClearObjective()
        obj = []
        for reac in reacs:
            back = self._BackReac(reac)
            indices = [reac]
            values = [1]
            if back != None:
                indices.append(back)
                values.append(-1)
            self.AddVar(obj=[1], types=[self.variables.type.binary], names=[reac+"_zero"])
            self.indicator_constraints.add(lin_expr=cplex.SparsePair(ind = indices, val = values), sense = 'E', rhs = 0.0, indvar = reac+"_zero", name = reac+"_isZero",complemented = 0)
#            self.indicator_constraints.add(lin_expr=cplex.SparsePair(ind = [reac], val = [1]), sense = 'G', rhs = 0.0, indvar = reac+"_pos", name = reac+"_is_pos",complemented = 0)
##            if back != None:
##                self.indicator_constraints.add(lin_expr=cplex.SparsePair(ind = [back], val = [1]), sense = 'E', rhs = 0.0, indvar = reac+"_neg", name = reac+"_not_neg",complemented = 1)
##                self.indicator_constraints.add(lin_expr=cplex.SparsePair(ind = [back], val = [1]), sense = 'G', rhs = 0.0, indvar = reac+"_neg", name = reac+"_is_neg", complemented = 0)
##            self.linear_constraints_names.append(reac+"_either_pos_or_neg")
##            self.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = [reac+"_neg",reac+"_pos"], val = [1,1])], senses=['L'], rhs=[1.0], names=[reac+"_either_pos_or_neg"])
#            obj.extend([reac+"])
#        self.SetObjective(obj)
        
    def createDiffObjective(self,reacs):
        """ reacs has to be a dictionary of reaction fluxes (w). all other fluxes are assumed to be 0"""
        self.ClearObjective()
        linpart = reacs.copy()
        #the direction is minimising!
        self.SetObjDirec("Min")
        for reac in reacs:
            back = self._BackReac(reac)
            if back != None:
                linpart[back] = linpart[reac]
            linpart[reac] = -linpart[reac]
            if not reac in self.variables_names:
                print reac + " not in the problem"
                continue
        #This means we get -w*x + x^2 as objective. which is, if minimized the same as minimizing the eukleadian norm.
        self.SetLinObjCoefsFromDic(linpart)
        self.SetQuadObjCoefsFromLists(self.variables_names,[0.5]*len(self.variables_names))
  
    def ClearObjective(self):
        self.objective.set_linear(zip(self.variables_names[:],[0]*len(self.variables_names)))
        self.objective.set_quadratic_coefficients(zip(self.variables_names[:],self.variables_names[:],[0]*len(self.variables_names)))
                                                  

    def SetObjCoefsFromDic(self, d, quadratic=False):
        if quadratic:
            self.SetQuadObjCoefsFromDic(d)
        else:
            self.SetLinObjCoefsFromDic(d)        

    def SetObjCoefsFromLists(self, cs, coefs, quadratic=False):
        if quadratic:
            self.SetQuadObjCoefsFromLists(cs, coefs)
        else:
            self.SetLinObjCoefsFromLists(cs, coefs)

    def SetObjDirec(self,direc="Min"):
        if not (direc=="Max" or direc=="Min"):
            raise exceptions.ValueError, direc
        self.objective.set_sense(ObjectiveMap[direc])

    def SetObjCoef(self,c,coef,quadratic=False):
        if quadratic:
            self.SetQuadObjCoef(c,coef)
        else:
            self.SetLinObjCoef(c,coef)

    def SetObjName(self,name):
        self.objective.set_name(name)

    def GetObjName(self):
        return self.objective.get_name()
    ##
    #   Access
    ##
    def GetObjective(self,Quadratic = False):
        if Quadratic:
            res = dict(zip(self.variables_names,self.objective.get_quadratic_coefficients(zip(self.variables_names,self.variables_names))))
            for r in res.keys():
                if res[r] == 0:
                    del res[r]
        else:
            res = dict(zip(self.variables_names,self.objective.get_linear()))
            for r in res.keys():
                if res[r] == 0:
                    del res[r]
        return res
    
    ##
    #   Quadratic
    ##
    
    def SetQuadObjective(self,reacs):
        
        if type(reacs)==types.DictType:
            target = {}
            for reac in reacs.keys():
                target[reac] = reacs[reac]
                back = self._BackReac(reac)
                if back != None:
                    target[back] = reacs[reac]
            self.objective.set_quadratic_coefficients(zip(target.keys(),target.keys(),target.values()))
            
        else:
            target = reacs[:]
            for reac in target:
                back = self._BackReac(reac)
                if back != None and back not in target:
                    target.append(back)
                
            self.objective.set_quadratic_coefficients(zip(target,target,[1]*len(target)))

    def SetQuadObjCoefsFromLists(self, cs, coefs):

        lenc = len(cs)
        if lenc != len(coefs):
            raise exceptions.IndexError, "cs and coefs not of same length !"

        self.objective.set_quadratic_coefficients(zip(cs,cs,coefs))

    def SetQuadObjCoefsFromDic(self, d):

        self.SetQuadObjCoefsFromLists(d.keys(), d.values())

    def SetQuadObjCoef(self, c, coef):
        """set a Quadratic coefficient for a specific variable"""
        if c in self.variables_names:
            self.objective.set_quadratic_coefficients(c,c,coef)
            
    ##
    #   Linear
    ##

    def SetLinObjective(self,reacs):
        #if we have a dictionary we 
        if type(reacs)==types.DictType:
            target = {}
            for reac in reacs.keys():
                target[reac] = reacs[reac]
                back = self._BackReac(reac)
                if back != None:
                    target[back] = reacs[reac]
            self.SetLinObjCoefsFromDic(target)

        else:
            target = reacs[:]
            temp = []
            for reac in target[:]:
                back = self._BackReac(reac)
                if not self.HasCol(reac):
                    raise exceptions.IndexError, "Column " + str(reac) + " not in the Problem"
                if back != None and back not in target:
                    target.append(back)
                else:
                    del back
                
            self.SetLinObjCoefsFromLists(target[:],[1]*len(target[:]))
            del target

    def SetLinObjCoefsFromLists(self, cs, coefs):
        lenc = len(cs)
        if lenc != len(coefs):
            raise exceptions.IndexError, "cs and coefs not of same length !"
        self.objective.set_linear(zip(cs[:], coefs[:]))


    def SetLinObjCoefsFromDic(self, d):
        self.SetLinObjCoefsFromLists(d.keys(), d.values())


    def SetLinObjCoef(self, c, coef):
        if c in self.variables_names:
            self.objective.set_linear(c,coef)
            
    
    ##
    ##
    #   Variable Modifications
    ##
    ##
        
    ##
    #   Variable Constraint Modifications
    ##
    def SetFixedFlux(self, reacs):

        for reac in reacs:

            rate = reacs[reac]
            if reac in self.sm.Reversed:
                rate *= -1
            self.SetColBounds(reac, rate, rate)
                    
    def SetColBounds(self, c, lo=None, hi=None):
        back = self._BackReac(c)
        if not c in self.variables_names:
            raise exceptions.ValueError, c + "is not part of this problem"
        if lo != None:
            if lo < 0 :
                if back == None:
                    raise exceptions.ValueError, "trying to set negative flux to irreversable reaction " + c
                else:
                    self.variables.set_upper_bounds(back,-lo)
                self.variables.set_lower_bounds(c,0)
            else:
                if back != None:
                    self.variables.set_upper_bounds(back,0)
                self.variables.set_lower_bounds(c,lo)
        else: 
            self.variables.set_lower_bounds(c,0)
            if back != None:
                self.variables.set_upper_bounds(back,self.infinity)
                
        if hi != None:
            if hi < 0 :
                if back == None:
                    raise exceptions.ValueError, "trying to set negative flux to irreversable reaction " + c
                else:
                    self.variables.set_lower_bounds(back,-hi)
                self.variables.set_upper_bounds(c,0)
            else:
                if back != None:
                    self.variables.set_lower_bounds(back,0)
                self.variables.set_upper_bounds(c,hi)
        else: 
            self.variables.set_upper_bounds(c,self.infinity)
            if back != None:
                self.variables.set_lower_bounds(back,0)

    def UnboundFlux(self, reac):

        self.SetColBounds(reac, 0, None)

        rback = self._BackReac(reac)
        if rback != None:
            self.SetColBounds(rback, 0, None)
                
    def FiniteBoundFlux(self, reac, lo, hi):

        if hi < lo:
            lo,hi = hi,lo
        self.SetColBounds(reac, lo,hi)

    def SetFluxBounds(self, reacs):
        """ Pre: reacs = {name:(lo,hi)...}, (lo<=hi OR hi==None) """

        for reac in reacs:

            lo, hi = reacs[reac]

            if lo == hi == None:
                self.UnboundFlux(reac)

            elif lo == hi:
                self.SetFixedFlux({reac:lo})

            elif not None in (lo,hi):
                self.FiniteBoundFlux(reac, lo, hi)

            else:                                  # lo == None xor hi == None
                rback = self._BackReac(reac)

                if rback == None:                  # irreversible
                    self.SetColBounds(reac, lo, hi)
                else:                              # reversible reaction
                    if lo == None:                 # lo == None
                        if hi < 0.0:
                            self.SetColBounds(reac, 0,0)
                            self.SetColBounds(rback, -hi, None)
                        else:
                            self.SetColBounds(reac,  0, hi)
                            self.SetColBounds(rback, 0, None)
                    else:
                        if lo < 0.0:                # hi == None
                            self.SetColBounds(reac, 0,None)
                            self.SetColBounds(rback, 0, -lo)
                        else:
                            self.SetColBounds(reac, lo,None)
                            self.SetColBounds(rback, 0, 0)
    ##
    #   Adding variables
    ##

    def AddCols(self,cols):

        tc = type(cols)
        if tc == int:
            if cols <1:
                raise exceptions.ValueError, "Can't add less then one col !"
            else:
                for i in range(cols):
                    name = "col_"+str(len(self.variables_names))
                    self.AddVar(names=[name])
#                    self.variables_names.append(name)

        elif tc == list:
            if True in map(self.HasCol,cols):
                raise exceptions.ValueError, "Can't add duplicate col index !"
            
            self.AddVar(names=cols)
#            self.variables_names.extend(cols)

        else:
            raise execeptions.TypeError, "cols must be int or list"

    def AddVar(self, obj=[], lb=[], ub=[], types='', names=[], columns=[]):
        if len(names) < 1:
            print "Need variable names"
            return
        else:
            for r in names:
                if r in self.variables_names:
                    raise exceptions.ValueError, "Duplicate variable name! " + r
                else:
                    self.variables_names.append(r)
                 
            self.variables.add(obj, lb, ub, types, names, columns)
            self.updateColDics()  # update self.cnames
            
    def AddIntegerVar(self,name,lower,upper):
        """Add an Integer Variable with the provided Name and upper and lower bounds to the problem"""
        self.AddVar(names=[name],lb=[lower],ub=[upper],types = ["I"])

    def AddBinaryVar(self,name):
        """Add a Binary Variable with the provided Name to the problem"""
        self.AddVar(names=[name],types = ["B"])
        
    ##
    #   Removing Variables
    ##
    
    def DelCols(self,cols):
        self.variables.delete(cols)
        for r in cols:
            del self.variables_names[self.variables_names.index(r)]
        self.updateColDics()  # update self.cnames
            
    ##
    #   Variable Modifications
    ##
    
    def SetColName(self,c,name):
        if HasCol(c):
            self.variables.set_names([(c,name)])
            self.variables_names[self.variables_names.index(c)] = name
            self.updateColDics()  # update self.cnames
        else:
            raise exceptions.IndexError, "Bad column index, " + str(c) + " not a variable name"

    ##
    ##
    #   Linear Constraint Modifications
    ##
    ##

    ##
    #   Adding Constraints
    ##
    def SetSumFluxConstraint(self,  reacs, total,  name):

        self.AddRows([name])
        reacs = reacs[:]
        for r in reacs[:]:
            rback = self._BackReac(r)
            if rback != None:
                reacs.append(rback)
        vals = [1]*len(reacs)
        self.linear_constraints.set_linear_components(name,[reacs,vals])
        self.SetRowBounds(name,  total, total)

        print "SetSumFluxCostraint testing !"

    def SetFluxConstraint(self,reaccoeffs,lb,ub,name):
        self.AddRows([name])
        reacs = reaccoeffs.keys()
        vals = []        
        for r in reacs[:]:
            rback = self._BackReac(r)
            if rback != None:
                    reacs.append(rback)
                    reaccoeffs[rback] = - reaccoeffs[r]
        for r in reacs:
            vals.append(reaccoeffs[r])
        self.linear_constraints.set_linear_components(name,[reacs,vals])
        self.SetRowBounds(name,lb,ub)
        

        
   
    def AddSumFluxConstraint(self,reacs,lo,hi,name):
        self.AddRows([name])
        reacs = reacs[:]
        for r in reacs[:]:
            rback = self._BackReac(r)
            if rback != None:
                reacs.append(rback)
        vals = [1]*len(reacs)
        self.linear_constraints.set_linear_components(name,[reacs,vals])
        self.SetRowBounds(name,  lo, hi)

    
    def AddRows(self,rownames):
        if True in map(self.HasRow,rownames):
            raise exceptions.ValueError, "Can't add duplicate row index !"
        self.AddLinConstraint(names=rownames)
#        self.linear_constraints.add(names=rownames[:])
#        self.linear_constraints_names.extend(rownames[:])
#        self.rnames[len(self.rnames)] = row   # add constraint in self.rnames


    def AddConstraint(self, name, lb , ub, lin_expr):
        """name is the name of the constraint,
           lb is the lower bound of the constraint
           up is the upper bound of the constraint
           lin_expr - is a list of tuples representing the constraint with each entry being: (var,coef)  e.g. [("Reaction1", 37)]
                      where var is the variable referenced and coef is the coefficient for the variable in this constraint.
        """
        if name =="":
            print "Need constraint name"
        else:
            if name in self.linear_constraints_names:
                    raise exceptions.ValueError, "Duplicate constraint name! " + r
            else:
                    flag,range,rhs = self.__getbnds(lb,ub)
                    ind,val = [],[]
                    for coef in lin_expr:
                        ind.append(coef[0])
                        val.append(coef[1])
                    #self.linear_constraints_names.append(r)
#                    print "name: " + str(len(name)) + "linexp" + str(len(lin_expr))
                    if flag == "FREE":
                        self.AddLinConstraint(lin_expr = [cplex.SparsePair(ind,val)], senses=['R'],  rhs = [cplex.infinity], range_values = [-2*cplex.infinity], names=[name])
                    else:
                        self.AddLinConstraint(lin_expr = [cplex.SparsePair(ind,val)], senses=[flag],  rhs = [rhs], range_values = [range], names=[name])
                    
    def AddFluxConstraint(self, name, lb , ub, lin_expr):
        """name is the name of the constraint,
           lb is the lower bound of the constraint
           up is the upper bound of the constraint
           lin_expr - is a list of tuples representing the constraint with each entry being: (var,coef)  e.g. [("Reaction1", 37)]
                      where var is the variable referenced and coef is the coefficient for the variable in this constraint.
        """
        if name =="":
            print "Need constraint name"
        else:
            if name in self.linear_constraints_names:
                    raise exceptions.ValueError, "Duplicate constraint name! " + r
            else:
                    flag,range,rhs = self.__getbnds(lb,ub)
                    ind,val = [],[]
                    for coef in lin_expr:
                        ind.append(coef[0])
                        val.append(coef[1])
                        back = self._BackReac(coef[0])
                        if back != None:
                            ind.append(back)
                            val.append(coef[1])
                            
                    #self.linear_constraints_names.append(r)
#                    print "name: " + str(len(name)) + "linexp" + str(len(lin_expr))
                    if flag == "FREE":
                        self.AddLinConstraint(lin_expr = [cplex.SparsePair(ind,val)], senses=['R'],  rhs = [cplex.infinity], range_values = [-2*cplex.infinity], names=[name])
                    else:
                        self.AddLinConstraint(lin_expr = [cplex.SparsePair(ind,val)], senses=[flag],  rhs = [rhs], range_values = [range], names=[name])
                           



    def AddLinConstraint(self, lin_expr=[], senses='', rhs=[], range_values=[], names=[]):
        if len(names) < 1:
            print "Need constraint names"
            return
        else:
            for r in names:
                if r in self.linear_constraints_names:
                    raise exceptions.ValueError, "Duplicate constraint name! " + r
                else:
                    self.linear_constraints_names.append(r)
            self.linear_constraints.add(lin_expr=lin_expr, senses=senses, rhs=rhs, range_values=range_values, names=names)
        self.updateRowDics() #update rnames
        
    ##
    #   Removing Constraints
    ##

    def DelRows(self,rows):
#        self.linear_constraints.delete(rows)
        for r in rows:
            self.DelRow(r)
#            if type(1) == type(r):
#                del self.linear_constraints_names[r]
#            if r in self.linear_constraints_names:
#                del self.linear_constraints_names[self.linear_constraints_names.index(r)]
#        updateRowDics() # update rnames

    def DelRow(self,id):
        self.linear_constraints.delete(id)
        if type(1) == type(id):
            del self.linear_constraints_names[id]
        if id in self.linear_constraints_names:
            del self.linear_constraints_names[self.linear_constraints_names.index(id)]
        self.updateRowDics()
                
    ##
    #   Modifing Constraints
    ##
    def SetRowVals(self,rowname,values):
        """pre: values is a list containing a number of values equal to the number of total variables in this QP."""
        if not rowname in self.linear_constraints_names:
            raise exceptions.IndexError, "Bad col index " + str(rowname) + " not a valid identifier"

        vals = []
        idx = []
        for i in range(len(values)):
            if values[i] != 0:
                vals.append(values[i])
                idx.append(i)
        self.linear_constraints.set_linear_components(rowname,[idx,vals])

    def SetRowFromFluxDic(self,rowname,reacs, lo,hi):
        vals = []
        idx = []
        rs = reacs.copy()
        for r in reacs:
            rback = self._BackReac(r)
            if rback != None:
                rs[rback] = reacs[r]
        if not rowname in self.linear_constraints_names:
            self.linear_constraints_names.append(rowname)
            self.linear_constraints.add(names=[rowname])
        self.linear_constraints.set_linear_components(rowname,[reacs.keys(),reacs.values()])
        self.SetRowBounds(rowname,  lo, hi)
        
    def SetRowBounds(self,rowname,lo=None,hi=None):
       #get upper and lower bounds
        flag,range,rhs = self.__getbnds(lo,hi)
        #if a constraint is to be freed it will be removed.
        if flag == "FREE":
            if rowname in self.linear_constraints_names:
                self.linear_constraints.delete(rowname)
                del self.linear_constraints_names[self.linear_constraints_names.index(rowname)]
        else:
            self.linear_constraints.set_rhs(rowname,rhs)
            self.linear_constraints.set_range_values(rowname,range)
            self.linear_constraints.set_senses(rowname,flag)

    def LoadMtxFromDic(self,dic, KeysAreRowNames=True):
        
        if KeysAreRowNames:
            
            for k in dic.keys():
                self.SetRowVals(k,dic[k])

    def LoadMtxFromLists(self, RowIdxs ,ColIdcs ,ElVals):

        #check if all rows and columns exist
        if not len(RowIdxs) == len(ColIdcs) == len(ElVals):
             raise exceptions.IndexError, "Inconsistant indices and/or values"

        for r in RowIdxs:
            if not r in self.linear_constraints_names:
                raise exceptions.IndexError, "Bad row index, " + str(r) + " not a constraint name"
        #remove all linear constraints
        self.linear_constraints.delete()
        self.linear_constraints_names = []
        
        for c in ColIdcs:
            if not c in self.variables_names:
                 raise exceptions.IndexError, "Bad column index, " + str(c) + " not a variable name"
        #set up the new constraints
        for r in RowIdxs:
            if not self.HasRow(r):
                self.linear_constraints.add(names=[r])
                self.linear_conatraints_names.append(r)
        
        for i in range(len(RowIdxs)):
            self.linear_constraints.set_linear_components(RowIdxs[i],[[ColIdcs[i]],[ElVals[i]]])

    def SetRowName(self,r,name):
        if HasRow(r):
            self.linear_constraints.set_names([(c,name)])
            self.linear_constraints_names[index(r)] = name
        else:
            raise exceptions.IndexError, "Bad row index, " + str(c) + " not a constraint name"
        self.updateRowDics()
    ##
    ##
    #   Quadratic constraint modifications
    ##
    ##

    ##
    #   Add quadratic constraints
    ##

    def AddQuadRow(self,quadvars,quadcoefs,name=None,lin_row=cplex.SparsePair(ind = [0], val = [0.0]),rhs=0,sense="L"):
        """Add a quadratic constraint,
            quadvars is a list of quadratic variables,
            quadcoefs is a list of the coefficients of the quadratic variables
            lin_row is a complete linear row with coefs for all variables
            this function  does not allow setting 'off diagonal' constraints"""
        if len(quadvars) != len(quadcoefs):
            raise exceptions.IndexError, "Inconsistant length of variables and coefficients"
        if name == None:
            name = "quad"+str(self.quadratic_constraints.get_num())
        if coefs != None:
            coefs = lin_row[:]
            lin_row=[[],[]]
            i = 0
            for coef in coefs:
                if coef != 0:
                    lin_row[0].append(i)
                    line_row[1].append(coef)
                i+=1       
        quad_expr=[[],[],[]]
        i = 0
        self.quadratic_constraints.add(name=name,lin_expr=lin_row,quad_expr=[quadvars,quadvars,quadcoefs],rhs=rhs,sense=sense)
        self.quadratic_constraints_names.extend(name)
        
    ##
    #   Delete Quadratic Constraints
    ##

    def DelQuadRow(self,r):
        if r in self.quadratic_constraints_names:
            self.quadratic_constraints.delete(r)
            del self.quadratic_constraints_names[self.quadratic_constraints_namesindex(r)]
        
        
    ##
    ##
    #   Solving and Solution access
    ##
    ##
    
    def Solve(self,PrintStatus=True):
        self.solve()
        if PrintStatus:
            print self.GetStatusMsg()    

    def GetStatusMsg(self):
        if self.solution.get_status() in Status2Str:    
            return Status2Str[self.solution.get_status()]
        else:
            return "CPLEX Solution Status Code: " + str(self.solution.get_status())
        
    def GetPrimSol(self, IncZeroes=False,FixBack=True,AsMtx=False,threshold=VerySmall, rounding = False):

        if self.GetStatusMsg()=="Optimal":
            sol = dict(zip(self.variables_names,self.solution.get_values()))
            cols = sol.keys()
            for r in cols:
                if IsZero(sol[r],threshold):
                    del sol[r]
        else:
            print self.GetStatusMsg()
            return {}
        for rev in Set.Intersect(self.sm.Reversed,  sol.keys()):
            sol[rev] *= -1

        if FixBack:
            for reac in sol.keys():
                if reac.endswith("_back"):
                    val = sol[reac]
                    if not IsZero(val,threshold):
                        #add up both values from revers and original reaction
                        r2 = reac[:-5]
                        if r2 in sol:
                            sol[r2]-= val
                        else:
                            sol[r2] = -val
                    del sol[reac]
        
        cols = sol.keys()
        for r in cols:
            if rounding:
                sol[r] = round(sol[r],int(round(-math.log(threshold,10))))
            if IsZero(sol[r],threshold):
                del sol[r]
                       
        if  AsMtx:
            rv = StoMat.StoMat(cnames=["lp_sol"],rnames=sol.keys(),Conv=float)
            for r in sol:
                rv[r][0] = sol[r]
        else:
            rv = sol

        return rv

    
    def GetRowDual(self,r):
        """get the Dual of constraint r """
        if self.GetStatusMsg() == "Optimal":
            return self.solution.get_dual_values(r)
        else:
            print "There is no optimal Solution"
            return None
    

    def GetColDual(self,c):
        """get the Dual of variable r IMPORTANT: This will not return the Flux value but the variable value (i.e. irrespectible of the reverse flux """
        if self.GetStatusMsg() == "Optimal":
            return self.solution.get_reduced_costs(c)
        else:
            print "There is no optimal Solution"
            return None
    

    def GetColPrimal(self, c):
        """get the Primal of variable r IMPORTANT: This will not return the Flux value but the variable value (i.e. irrespectible of the reverse flux """
        if self.GetStatusMsg() == "Optimal":
            return self.solution.get_values(c)
        else:
            print "There is no optimal Solution"
            return None
    

    def GetRowPrimal(self, r):
        """get the Primal of constraint r """
        if self.GetStatusMsg() == "Optimal":
            return self.solution.get_linear_slacks(r)
        else:
            print "There is no optimal Solution"
            return None
    

    def GetObjDir(self):
        """get the Directionality of the objective (either Min or Max)"""
        return ObjectMap[self.objective.get_sense()]

    def GetObjVal(self):
        """get the Objective Value if there is a valid solution of constraint r """
        if self.GetStatusMsg() == "Optimal":
            return self.solution.get_objective_value()
        else:
            print "There is no optimal Solution"
            return None


    def GetFluxValue(self,flux):
        if self.GetStatusMsg() == "Optimal":
            val = self.solution.get_values(flux)
            back = _BackReac(flux)
            if back != None:
                val -= self.solution.get_values(back)
            return val
        else:
            print "There is no optimal Solution"
            return None
       
    
    ##
    ##
    #   General Access and Information Methods
    ##
    ##

    ##
    #   Matrix Information Methods
    ##
    def GetSimplexAsMtx(self):

        rnames = self.GetRowNames()
        cnames = self.GetColNames()
        rows   = self.GetRowsAsLists()

        rv = DynMatrix.matrix(rnames=rnames, cnames=cnames, Conv=float)
        rv.rows =rows

        return rv
    def GetRowsAsLists(self):
        #create an empty matrix
        matrix = []
        for row in self.linear_constraints_names:
            matrix.append(self.GetRow(row)[:])
        return matrix

    ##
    #   Variable/Constraint info retrieving
    ##
    def GetColNames(self):

        return self.variables_names

    def GetColIdx(self,c):
        return self.variables_names.index(c)
    
    def GetColIdxs(self,cols):
        return map(self.variables_names.index,cols) # incompatible with ScrumPyLP
    

    def GetCol(self,c):
        col = []
        cidx = self.GetColIdx(c)
        for r in self.linear_constraints.get_rows():
            if cidx in r.ind:
                col.append(r.val[r.ind.index(cidx)])
            else:
                col.append(0)
                
        return col
    
    def GetRowNames(self):
        return self.linear_constraints_names

    def GetRowIdx(self,r):
        return self.linear_constraints_names.index(r)
    
    def GetRowIdxs(self,rows):
        return map(self.linear_constraints_names.index,rows) # incompatible with ScrumPyLP
 
    def GetRow(self,row):
        comps = self.linear_constraints.get_rows()[self.linear_constraints.get_indices(row)]
        rv = [0.0]*len(self.variables_names)
        for j in range(len(comps.ind)):
            rv[comps.ind[j]] = comps.val[j]
            
        return rv
    
    ##
    #   General property retrieving
    ##
    
    
    def GetNumCols(self):
        return len(self.variables_names)

    def GetNumRows(self):
        return len(self.linear_constraints_names)

    def HasRow(self,rowname):
        if rowname in self.linear_constraints_names:
            return True
        else:
            return False

    def HasCol(self, c):
        """pre: True
          post: HasCol(c) => c is a valid col index """
        if c in self.variables_names:            
            return True
        else:
            return False

    def GetColName(self, c):
        #must be smaller then the number of variables, larger than zero and a integer)
        if c < self.variables.get_num() and c >= 0 and type(c) == type(1):
            return self.variables_names[c]

    def GetRowName(self, r):
        #must be smaller then the number of variables, larger than zero and a integer)
        if r < self.linear_constraints.get_num() and r >= 0 and type(c) == type(1):
            return self.linear_constraints_names[c]
        
    def GetDims(self):
        return (self.linear_constraints.get_num(), self.variables.get_num())

    def SetName(self,name):
        self.set_problem_name(name)

    def GetName(self):
        return self.get_problem_name()
                
    def PrintMtx(self):
        """ print current constraint matrix on stdout """

        rv =""
        nr,nc = self.GetDims()
        print " "*10,
        for c in range(nc):
            print self.GetColName(c).ljust(6),
        print
        for r in range(nr):
            print self.GetRowName(r).ljust(10), " ".join(map(lambda x:("%3.3g"%x).ljust(6), self.GetRow(r)))
            
    ##
    ##
    #   Helper Functions
    ##
    ##
    def __getbnds(self,lo,hi):
        if   lo == None and hi == None:
            flag, range, rhs = "FREE", 0.0, 0.0
        elif lo != None and hi == None:
            flag, range, rhs = GREATER_THAN, 0.0, lo
        elif lo == None and hi != None:
            flag, range, rhs = LESS_THAN, 0.0, hi
        elif lo != None and hi != None:
            if lo < hi:
                flag, range, rhs = RANGED, -(hi-lo) , hi
            elif lo==hi:
                flag, range, rhs = EQUAL_TO, 0.0, hi
            else:
                raise exceptions.ValueError, "lo > hi !"

        return flag, range, rhs


        
    def _BackReac(self,reac):
        if reac in self.sm.Backs:
            return self.sm.Backs[reac]
        return None
 
   
    def FixedFluxScan(self, reaction, lo, hi, n_p):

        rv = DataSets.DataSet(ItemNames=[reaction,"ObjVal"])
        lo = float(lo)
        hi = float(hi)
        inc = (hi - lo)/(n_p-1)
        cur = lo

        for n in range(n_p):
            self.SetFixedFlux({reaction:cur})
            self.Solve(False)

            sol = self.GetPrimSol()

            if len(sol)==0:
                obval = float("NaN") # indicate fail
            else:
                obval =  self.GetObjVal()

            sol["ObjVal"] = obval
            sol[reaction] = cur
            rv.UpdateFromDic(sol)

            cur += inc

        return rv

    ##
    # Add ATP Consumption and production variables
    ##
    
    def AddATP(self,singleprod,doubleprod,singlecons,doublecons,revs):
        """ singleprod are ones which produce a single ATP from ADP or ADP from AMP
        doubleprod are reactions adding diphosphate"""
        atpproddic = {}
        atpcondic = {}
        for r in singleprod:
            atpproddic[r] = 1
            if r in revs:   
                atpcondic[self._BackReac(r)] = 1
        for r in doubleprod:
            if r in atpproddic:
                atpproddic[r] +=2
            else:
                atpproddic[r] = 2
            if r in revs:
                back = self.BackReac(r)
                if back in atpcondic:
                    atpcondic[back] +=2
                else:
                    atpcondic[back] = 2
        for r in singlecons:
            if r in atpcondic:
                atpcondic[r] +=1
            else:
                atpcondic[r] = 1
            if r in revs:
                back = self.BackReac(r)
                if back in atpproddic:
                    atpproddic[back] +=1
                else:
                    atpproddic[back] = 1
            
        for r in doublecons:
            if r in atpcondic:
                atpcondic[r] +=2
            else:
                atpcondic[r] = 2
            if r in revs:
                back = self.BackReac(r)
                if back in atpproddic:
                    atpproddic[back] +=2
                else:
                    atpproddic[back] = 2
           
    
