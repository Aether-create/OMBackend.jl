#= The ModelingToolkit.jl package is licensed under the MIT "Expat" License:
# Copyright (c) 2018-20: Christopher Rackauckas, Julia Computing.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  SOFTWARE.
=#


#=
* This file is part of OpenModelica.
*
* Copyright (c) 2021-CurrentYear, Open Source Modelica Consortium (OSMC),
* c/o Linköpings universitet, Department of Computer and Information Science,
* SE-58183 Linköping, Sweden.
*
* All rights reserved.
*
* THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
* THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
* ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
* RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
* ACCORDING TO RECIPIENTS CHOICE.
*
* The OpenModelica software and the Open Source Modelica
* Consortium (OSMC) Public License (OSMC-PL) are obtained
* from OSMC, either from the above address,
* from the URLs: http:www.ida.liu.se/projects/OpenModelica or
* http:www.openmodelica.org, and in the OpenModelica distribution.
* GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
*
* This program is distributed WITHOUT ANY WARRANTY; without
* even the implied warranty of  MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
* IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
*
* See the full OSMC Public License conditions for more details.
*
=#


#=
This file contains "hacks".
This is done in order
to get the equations on a MTK compatible format before calling functions such as structurally simplify.
TODO:
!Adjust the uncessary string conversions!
=#

#=
  So we know about t an der in the global scope.
  This is needed for the rules below to match correctly.
=#
ModelingToolkit.@variables t
const D = Differential(t)

"""
  Temporary rewrite function. Not very pretty...
  Original code by Chris R. Expanded to fix terms of type X * D(Y).
  The solution to solve it is not pretty and is probably flaky.
#istree returns true if x is a term. If true, operation, arguments must also be defined for x appropriately.
"""
function move_diffs(eq::Equation; rewrite)
  # Do not modify `D(x) ~ ...`, already correct
  res =
    if !(istree(eq.lhs) && operation(eq.lhs) isa Differential) && !(istree(eq.lhs) && operation(eq.lhs) isa Difference)
      local _eq
      _eq = eq.rhs-eq.lhs
      rhs = rewrite(_eq)
      if rhs === nothing
        return eq
      end
      lhs = _eq - rhs
      if !(lhs isa Number) && (operation(lhs) isa Differential)
        lhs ~ -rhs
      elseif !(lhs isa Number) && (operation(lhs) == *)
        #=
          This code is probably quite flaky, however, it should not be needed after similar things are introduced in MTK.
          TODO: Handle this more elegantly using rewrite rules instead.
        =#
        local newRhs
        local newLhs
        for arg in arguments(lhs)
          local isTermAndNotDifferential = istree(arg) && !(operation(arg) isa Differential)
          local argIsANumberOrSymbolButNotTerm = (arg isa Number || (arg isa SymbolicUtils.BasicSymbolic{Real} && !istree(arg)))
          if argIsANumberOrSymbolButNotTerm || isTermAndNotDifferential
            newRhs = substitute(rhs, rhs => rhs / arg)
            rhs = newRhs
            newLhs = substitute(lhs, lhs => lhs / arg)
            lhs = newLhs
            #@info "New equation in the for loop is $(lhs) = $(rhs)"
          end
        end
        tmp = ~(newLhs, -newRhs)
        tmp
      else
        -lhs ~ rhs
      end
    else
      eq
    end
  return res
end

"""
  Rewrite equations that do not conform to the requirements of MTK
"""
function rewriteEquations(edeqs, iv, eVars, ePars, simCode)
  # println("Recived #edeqs")
  # println(length(edeqs))
  local der = ModelingToolkit.Differential(t)
  #= Remove the t's =#
  eVars = [Symbol(replace(string(i), "(t)" => "")) for i in eVars]
  #Temporary fix for the ESCIMO climate model: eVars = vcat(eVars, [Symbol("combi_Population_Lookup_bn_y")])
  preEval = quote
    vars = ModelingToolkit.@variables begin
      $(eVars...)
    end
    pars = ModelingToolkit.@parameters begin
      $(ePars...)
    end
  end
  #Hardcoded
  #= Make the derivative symbol known =#
  eval(preEval)
  eval(:(der = ModelingToolkit.Differential(t)))
  eval(:(import ModelingToolkit.IfElse))
  #= Make the external runtime available if it is used. =#
  if simCode.externalRuntime
    eval(generateExternalRuntimeImport(simCode))
  end
  #=
  Register all functions s.t they can be used by the symbolic transformations.
  =#
  #= Delay evaluation of the register expression until we have constructed the functions =#
  local ftrs = generateRegisterCallsForCallExprs(simCode; funcArgGen = AlgorithmicCodeGeneration.generateIOL)
  for f in ftrs
    eval(f)
  end
  #==#
  local deqs = evalEDeqs(edeqs)
  #println(debugRewrite(deqs, iv, vars, parameters; separator="\n"))
  #= Rewrite equations =#
  D = Differential(iv)
  local r1 = SymbolicUtils.@rule ~~a * D(~~b) * ~~c => 0
  local r2 = SymbolicUtils.@rule D(~~b) => 0
  local remove_diffs = SymbolicUtils.Postwalk(SymbolicUtils.Chain([r1,r2]))
  local usedStates = Set()
  local rewrittenDeqs = Symbolics.Equation[]
  local req
  for eq in deqs
    #= Only do the rewrite for the differentials. The others have already been rewritten.=#
    eqStr = string(eq)
    if (contains(eqStr, "Differential")) # TODO expensive comp! Needs to be optimized.
    req = move_diffs(eq, rewrite = remove_diffs)
      #@info "Left hand side of the equation" req.lhs
      if req.lhs isa Real
        push!(rewrittenDeqs, req)
      elseif !(req.lhs in usedStates)
        #@info "Not a duplicate" req.lhs
        #@info "Used equations are" usedStates
        push!(rewrittenDeqs, req)
        push!(usedStates, req.lhs)
      else
        #@info "Duplicate:" req.lhs
        push!(rewrittenDeqs, eq)
      end
    else
      push!(rewrittenDeqs, eq)
    end
  end
  #OMBackend.debugWrite("codeAfterBackendPreprocessing.log", debugRewrite(rewrittenDeqs, iv, vars, parameters; separator="\n"))
  return rewrittenDeqs
end

"""
This function evaluates the supplied equations.
In the case we are unable to evaluate them, we currently hack it by some string conversions.
TODO:
  Fix me do this the proper way.
  This routine is way way to slow currently...
"""
function evalEDeqs(edeqs)
  local deqs = []
  for e in edeqs
    try
      eq = eval(e)
      if typeof(eq.lhs) == Int64 && eq.lhs == 0
        push!(deqs, eq)
      else
        local unSimplifiedString = string(e)
        unSimplifiedString = replace(unSimplifiedString, "&&" => "&")
        unSimplifiedString = replace(unSimplifiedString, r"\bbegin\b" => "(")
        unSimplifiedString = replace(unSimplifiedString, r"\bend\b" => ")")
        #println(unSimplifiedString)
        estrExp = Meta.parse(unSimplifiedString)
        estrExp2 = eval(estrExp)
        estrExp2LHS = estrExp2.lhs
        @assign estrExp2.lhs = 0
        @assign estrExp2.rhs = estrExp2.rhs - estrExp2LHS
        push!(deqs, estrExp2)
      end
    catch ex
      #Hack
      local unSimplifiedString = string(e)
      unSimplifiedString = replace(unSimplifiedString, "&&" => "&")
      unSimplifiedString = replace(unSimplifiedString, r"\bbegin\b" => "(")
      unSimplifiedString = replace(unSimplifiedString, r"\bend\b" => ")")
      #println(unSimplifiedString)
      estrExp = Meta.parse(unSimplifiedString)
      estrExp2 = eval(estrExp)
      estrExp2LHS = estrExp2.lhs
      @assign estrExp2.lhs = 0
      @assign estrExp2.rhs = estrExp2.rhs - estrExp2LHS
      push!(deqs, estrExp2)
    end
  end
  return deqs
end

function debugRewrite(deqs, t, vars, parameters; separator = ",")
  local buffer = IOBuffer()
  print(buffer, "@variables t;")
  print(buffer, "vars2 = @variables")
  print(buffer, "(")
  for v in vars
    print(buffer, v, separator)
  end
  print(buffer, ");")
  print(buffer, "eqs2 =")
  print(buffer, "[")
  for eq in deqs
    print(buffer, replace(string(eq), "(t)" => "", "Differential" => "D"), separator)
  end
  print(buffer, "];")
  return String(take!(buffer))
end

rewriteEq(eq) = begin
  local eqStr = string(eq)
  res = Meta.parse(replace(eqStr, "Differential(t)" => "D"))
  res
end

"""
$(SIGNATURES)

Structurally simplify algebraic equations in a system and compute the
topological sort of the observed equations. When `simplify=true`, the `simplify`
function will be applied during the tearing process. It also takes kwargs
`allow_symbolic=false` and `allow_parameter=true` which limits the coefficient
types during tearing.

The optional argument `io` may take a tuple `(inputs, outputs)`.
This will convert all `inputs` to parameters and allow them to be unconnected, i.e.,
simplification will allow models where `n_states = n_equations - n_inputs`.
"""
# function structural_simplify(sys::ModelingToolkit.AbstractSystem, io = nothing; simplify = false, kwargs...)
#   @info "Calling custom structural_simplify"
#   sys = expand_connections(sys)
#   state = TearingState(sys)
#   has_io = io !== nothing
#   has_io && markio!(state, io...)
#   state, input_idxs = ModelingToolkit.inputs_to_parameters!(state, io)
#   sys, ag = ModelingToolkit.alias_elimination!(state; kwargs...)
#   check_consistency(state, ag)
#   sys = dummy_derivative(sys, state, ag; simplify)
#   fullstates = [map(eq -> eq.lhs, observed(sys)); states(sys)]
#   @set! sys.observed = ModelingToolkit.topsort_equations(observed(sys), fullstates)
#   ModelingToolkit.invalidate_cache!(sys)
#   return has_io ? (sys, input_idxs) : sys
# end



# function structural_simplify(sys::ModelingToolkit.AbstractSystem, io = nothing; simplify = false, kwargs...)
#   @info "Calling custom structural_simplify"
#   #sys = ModelingToolkit.ode_order_lowering(sys)
#   #sys = ModelingToolkit.dae_index_lowering(sys)
#   #sys = ModelingToolkit.tearing(sys; simplify = simplify)
#   sys = ModelingToolkit.structural_simplify(sys, simplify = simplify)
#   return sys
# end

"""
TODO:
Document why some parts here are outcommented
The irreductable variables scheme does not work using plain simplify.


It should be noted that for some models both running tearing and structurally simplify is needed.
Report and issue for the MTK reporters giving an example of this behavior.

One example is running tearing twice broke the system

"""
function structural_simplify(sys::ModelingToolkit.AbstractSystem,
                             io = nothing;
                             simplify = false,
                             kwargs...)
  #@info "Calling custom structural_simplify"
  #sys = ModelingToolkit.ode_order_lowering(sys)
  #sys = ModelingToolkit.dae_index_lowering(sys)
  #sys = ModelingToolkit.tearing(sys; simplify = simplify)
   if simplify #Note report this to the developers of modeling toolkit.
     sys = ModelingToolkit.ode_order_lowering(sys)
     sys = ModelingToolkit.dae_index_lowering(sys)
     #sys = ModelingToolkit.tearing(sys; simplify = false)
     #Note some system breaks if tearing is run twice.
     # Note2 In some cases we need to do index reduction before simplify
   end
  sys = ModelingToolkit.structural_simplify(sys, simplify = false)
  return sys
 end

function getSyms(odeFunc::ODEFunction)
  return odeFunc.syms
end
