module RuntimeUtil

import Absyn
import ListUtil
import ModelingToolkit
import OMBackend
import OMBackend.SimulationCode
import OMFrontend
import OMFrontend.Main
import OMFrontend.Main.SCodeUtil
import OMFrontend.Main.Util
import SCode

using MetaModelica

"""
  Wrapper to a function in SCode util.
  inIdent is a string since frontend automatically remove certain parameters.
  TODO: Also search for components in the innermost class.
"""
function getElementFromSCodeProgram(inIdent::String, inClass::SCode.Element)
  result = SCodeUtil.getElementNamed(inIdent, inClass)
  return result
end

"""
  Given a name sets that element to a new value.
  It then returns the modified SCodeProgram.
Currently254 it is assumed to be at the top level of the class.
TODO: Fix for sublevels as well
"""
function setElementInSCodeProgram!(inIdent::String, newValue::T, inClass::SCode.Element) where {T}
  #= Get all elements from the class together with the corresponding names =#
  local elementToReplace = getElementFromSCodeProgram(inIdent, inClass)::SCode.Element
  local elements = listArray(SCodeUtil.getClassElements(inClass))::Vector{SCode.Element}
  local i = 1
  local indexOfElementToReplace = 0
  local outClass = inClass
  for element in elements
    if SCodeUtil.elementNameEqual(element, elementToReplace)
      indexOfElementToReplace = i
      break
    end
    i += 1
  end
  local modification = SCodeUtil.getComponentMod(elementToReplace)
  @assign modification.binding = makeCondition(newValue)
  @assign elementToReplace.modifications = modification
  elements[indexOfElementToReplace] = elementToReplace
  @assign outClass.classDef.elementLst = arrayList(elements)
  return outClass
end

makeCondition(val::Bool) = begin
  SOME(Absyn.BOOL(val))
end
makeCondition(val::Int) = begin
  SOME(Absyn.INTEGER(val))
end
makeCondition(val::Real) = begin
  SOME(Absyn.REAL(string(val)))
end
makeCondition(val) = begin
  throw("Only primitive values {Integer, Boolean, Real} are currently supported in a Recompilation call")
end

"""
  Converts a symbol (of a MTK variable) to a string.
"""
function convertSymbolToString(symbol::Symbol)
  res = replace(String(symbol), "(t)" => "")
  #= Remove prefixes in front of variables =#
  return res
end

"""
  Converts a list of symbols to a list of strings
"""
function convertSymbolsToStrings(symbols::Vector{Symbol})
  map(convertSymbolToString, symbols)
end

"""
  This function maps variables between two models during a structural change with recompilation.
  It returns a new vector of u₀ variables to initialize the new model.
  We do so by assigning the old values when the structural change occured for all variables
  that occured in the model before the structural change.
TODO:
Remove the special case.
The reason for it is in some example the containing model changes name wheras in others it does not.
Furthermore, the way indices are handled need to be fixed since the index after reduction in MTK is not the same as
the index in the simulation code stage of the backend.
For now we return the array for the special case with dynamic overconstrained connectors.
"""
function createNewU0(symsOfOldProblem::Vector{Symbol},
                     symsOfNewProblem::Vector{Symbol},
                     newHT,
                     initialValues,
                     integrator,
                     specialCase)
#  @debug "Length of old and new" length(symsOfOldProblem) length(symsOfNewProblem)
  local newU0 = Float64[last(initialValues[idx]) for idx in 1:length(symsOfNewProblem)]
  local variableNamesOldProblem = RuntimeUtil.convertSymbolsToStrings(symsOfOldProblem)
  local variableNamesNewProblem = RuntimeUtil.convertSymbolsToStrings(symsOfNewProblem)
  @debug variableNamesOldProblem
  @debug variableNamesNewProblem
  #= It was assumed to only be real variable not discretes, which might have other indices? =#
  local variableNamesWithoutPrefixesOP
  local variableNamesWithoutPrefixesNP
  if ! specialCase
    variableNamesWithoutPrefixesOP = [replace(k, r".*_" => "") for k in variableNamesOldProblem]
    variableNamesWithoutPrefixesNP = [replace(k, r".*_" => "") for k in variableNamesNewProblem]
  else
    return integrator.u
  end
  local largestProblem = if length(variableNamesOldProblem) > length(variableNamesOldProblem)
    variableNamesWithoutPrefixesOP
  else
    variableNamesWithoutPrefixesNP
  end
  for v in largestProblem
    local varNameWithoutPrefix = v
    if varNameWithoutPrefix in variableNamesWithoutPrefixesOP && varNameWithoutPrefix in variableNamesWithoutPrefixesNP
      local oldIndices = findall((x)-> x == varNameWithoutPrefix, variableNamesWithoutPrefixesOP)
      @assert(length(oldIndices) == 1, "Zero or more than one variable with that name. Size of oldIndicies was $(length(oldIndices)). Name was $(v)")
      idxOldVar = first(oldIndices)
      #= Locate the index of a variable with that name in the set of new variables=#
      local indices = findall((x)-> x == varNameWithoutPrefix, variableNamesWithoutPrefixesNP)
      #= I assume here that there are no duplicate variables =#
      @assert(length(indices) == 1, "Zero or more than one variable with that name. Size of indices was $(length(indices)). Name was $(v)")
      local varNameNew = variableNamesNewProblem[first(indices)]
      local idxNewVar = getIdxFromEntry(newHT[varNameNew])
      newU0[idxNewVar] = integrator.u[idxOldVar]
    end
  end
  return newU0
end

"""
  Gets the index from an entry to the symbol table
"""
function getIdxFromEntry(entry::Tuple)::Int
  first(entry)
end


"""
  Creates a new flat model.
  This model either has the equation specified in the if or these equations are to be added to the model
  when this condition is true.
  The name of this model
"""
function createNewFlatModel(flatModel,
                            unresolvedEquations,
                            newEquations)
  local newFlatModel =
    OMFrontend.Main.FLAT_MODEL(flatModel.name,
                               flatModel.variables,
                               flatModel.equations,
                               flatModel.initialEquations,
                               flatModel.algorithms,
                               flatModel.initialAlgorithms,
                               nil,
                               NONE(),
                               nil,
                               nil,
                               Bool[], #= TODO: This equation might need to be changed. =#
                               flatModel.comment)
  println("Adding equations")
  println("********************************************************************")
  @debug length(newEquations)
  for e in newEquations
    println(OMFrontend.Main.toString(e))
  end
  println("********************************************************************")

  # println("Existing equations:")
  # println("********************************************************************")
  # @debug "Length of existing system:" length(newFlatModel.equations)
  # for e in newFlatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  @assign newFlatModel.equations = listAppend(unresolvedEquations, newEquations)
  # println("Unresolved System:")
  # println("********************************************************************")
  # for e in newFlatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  #=
    1. Reresolve the connect equations.
    2. Perform constant evaluation.
    3. Run the simplification pass.
  =#
  newFlatModel = OMFrontend.Main.resolveConnections(newFlatModel, newFlatModel.name)
  newFlatModel = OMFrontend.Main.evaluate(newFlatModel)
  newFlatModel = OMFrontend.Main.simplifyFlatModel(newFlatModel)
  # println("Final System")
  # @debug "Length of final system:" length(newFlatModel.equations)
  # println("********************************************************************")
  # for e in newFlatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  return newFlatModel
end

"""
  Creates a new flat model with a set of connection equation removed.
  Note that the flat model passed to this function does not have any active equations.
  This means that this flat model is a new flat model with the unbreakable branches removed.
"""
function createNewFlatModel(flatModel,
                            idx::Int,
                            unresolvedEquations)
  local aDoccs = flatModel.active_DOCC_Equations
  aDoccs[idx] = false
  local newFlatModel =
    OMFrontend.Main.FLAT_MODEL(flatModel.name,
                               flatModel.variables,
                               flatModel.unresolvedConnectEquations,
                               flatModel.initialEquations,
                               flatModel.algorithms,
                               flatModel.initialAlgorithms,
                               nil,
                               NONE(),
                               flatModel.DOCC_equations,
                               flatModel.unresolvedConnectEquations,
                               aDoccs,
                               flatModel.comment)
  println("Create new model")
  local variablestoReset = resolveDOOCConnections(flatModel, flatModel.name)
  println("********************************************************************")
  println("Existing equations:")
  println("********************************************************************")
  @info "Length of OLD FLAT MODEL:" length(OMBackend.CodeGeneration.OLD_FLAT_MODEL.equations)
  for e in OMBackend.CodeGeneration.OLD_FLAT_MODEL.equations
    println(OMFrontend.Main.toString(e))
  end

  # println("********************************************************************")
  # println("New System:")
  # println("********************************************************************")
  # @info "Length of new system:" length(newFlatModel.unresolvedConnectEquations)
  # for e in newFlatModel.unresolvedConnectEquations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  #=
    1. Reresolve the connect equations.
    2. Perform constant evaluation.
    3. Run the simplification pass.
  =#
  newFlatModel = OMFrontend.Main.resolveConnections(newFlatModel, newFlatModel.name)
  newFlatModel = OMFrontend.Main.evaluate(newFlatModel)
  newFlatModel = OMFrontend.Main.simplifyFlatModel(newFlatModel)
  println("Final System")
  @info "Length of final system:" length(newFlatModel.equations)
  println("********************************************************************")
  for e in newFlatModel.equations
    println(OMFrontend.Main.toString(e))
  end
  println("********************************************************************")
  return newFlatModel
end

"""
  Resolves the system at the time of the structural change.
"""
function resolveDOOCConnections(flatModel, name)
  #= Get the relevant OCC graph =#
  local (searchGraph, rootVariables, rootEquations) = SimulationCode.getOCCGraph(flatModel)
  local pathsForRoots = Dict{String, Vector{String}}()
  for rv in rootVariables
    p = findPath(searchGraph, rv)
    pathsForRoots[OMFrontend.Main.toString(rv)] = p
  end
  for key in keys(pathsForRoots)
    println("Assignments for $(key):")
    vars = pathsForRoots[key]
    for v in vars
      println("$(v) := $(key)")
    end
  end
  local rootSources = Dict{String, String}()
  #= These are the equations for which the chain starts =#
  for (lhs, rhs) in rootEquations
    rootSources[OMFrontend.Main.toString(lhs)] = OMFrontend.Main.toString(rhs)
  end
  return (pathsForRoots, rootSources)
end

"""
author:johti17
Iterative DFS:
  Finds the path for a root variable passed as inV
"""
function findPath(g::Dict{String, Vector{String}}, inV)
  local v = OMFrontend.Main.toString(inV)
  local S = String[]
  local discovered = String[] #= Should ideally be int instead... =#
  push!(S, v)
  while  ! isempty(S)
    local v = pop!(S)
    if ! (v in discovered)
      push!(discovered, v)
      neighbours = g[v]
      for n in neighbours
        push!(S, n)
      end
    end
  end
  return discovered[2:end]
end

"""
Temporary function.
Evaluates a discrete events
Assuming one variable and one event that is changed.
TODO: Generalize later
"""
function evalDiscreteEvents(discreteEvents, u, t, system)
  local events = Tuple{Int, Bool, Bool}[]
  for de in discreteEvents
    push!(events, evalDiscreteEvent(de, u, t, system))
  end
  events = filter((x) -> last(x), events)
  return events
end

"""
TODO:
Refactor this function
"""
function evalDiscreteEvent(discreteEvent, u, time, system)
  @assert(length(discreteEvent.affects ) == 1, "Only length one of discrete affects supported")
  @info "New iteration\n\n\n\n"
  local affect = first(discreteEvent.affects)
  local condition = discreteEvent.condition
  local args = condition.arguments
  local operator = condition.f
  local stateVars = ModelingToolkit.states(system)
  local lhs = first(args)
  local rhs = last(args)
  local lhsIdx::Int = 0
  local rhsIdx::Int = 0
  local isChanged = false
  local varDeps = Int[]
  #Assuming a ! for this case
  local shouldApplyNegation = false
  if length(args) == 1
    lhs = first(first(args).arguments)
    rhs = last(first(args).arguments)
    shouldApplyNegation = true
    operator = first(args).f
  end
  if typeof(lhs) != Float64 && string(lhs) != string(system.iv)
    @info "args" args
    @info "lhs" lhs
    @info "lhs" typeof(lhs)
    @info "rhs" rhs
    lhsIdx = findfirst((x)->x==1, indexin(stateVars, [lhs]))
  end
  if typeof(rhs) != Float64 && string(rhs) != string(system.iv)
    rhsIdx = findfirst((x)->x==1, indexin(stateVars, [rhs]))
  end
  rhsValue = if rhsIdx != 0
    varDeps = getVariableEqDepedenceViaIdx(rhsIdx, system)
    @info "varDeps rhs" varDeps
    rootIdx = getRootEquation(varDeps)
    getConstantValueOfEq(rootIdx, system)
  elseif string(rhs) == string(system.iv)
    time
  else
    rhs
  end
  lhsValue = if lhsIdx != 0
    varDeps = getVariableEqDepedenceViaIdx(lhsIdx, system)
    @info "varDeps lhs" varDeps
    @info "root eq" getRootEquation(varDeps)
    rootIdx = getRootEquation(varDeps)
    getConstantValueOfEq(rootIdx, system)
  elseif string(lhs) == string(system.iv)
    time
  else
    lhs
  end
  local affectIdx = findfirst((x)->x==1, indexin(stateVars, [affect.lhs]))
  local affectNewValue = affect.rhs
  @info "lhs value was" lhsValue
  @info "rhs value was" rhsValue
  if shouldApplyNegation
    if operator(lhsValue, rhsValue) == false
      @info operator(lhsValue, rhsValue)
      #= Also assuming here that the lhs is a variable and the rhs is a value =#
      @info "Branch 1 Value was changed" discreteEvent.condition.f
      isChanged = true
    end
  else
    if operator(lhsValue, rhsValue)
      isChanged = true
    end
  end
  return (affectIdx, affect.rhs, isChanged)
end

"""
 Given a variable index, returns the equations that the variable at this index is dependent on.
  That is, equations in which this variable is referenced.
"""
function getVariableEqDepedenceViaIdx(idx::Int, system)
  #= Get all equation dependencies for the current system =#
  local equationDependencies = ModelingToolkit.equation_dependencies(OMBackend.Runtime.REDUCED_SYSTEM)
  local vars = ModelingToolkit.states(OMBackend.Runtime.REDUCED_SYSTEM)
  local totalDependencies = Int[]
  #= Go through each equation =#
  for (equationIndex, equationDep) in enumerate(equationDependencies)
    #= Skip equations without dependencies =#
    if isempty(equationDep)
      continue
    end
    #=
      If the equation dependency is not empty
      Check if it depends on our variable
    =#
    if first(indexin([vars[idx]], equationDep)) !== nothing
      #=
      In this case we know that this equation is a depedency of the supplied variable.
      Add this equation as a possible depdency to totalDependencies
      =#
      push!(totalDependencies, equationIndex)
    end
  end
  #= We now have All indices of the variables our equation depends on =#
  return totalDependencies
end

"""
  Get the top level equation if such equation exist for a given set of equations
  Note that if the supplied equation is not solved at the top level, this function returns 0
"""
function getRootEquation(equationIndices; usedEqIndices = Set())::Int
  @info "New call" equationIndices
  local G = ModelingToolkit.asgraph(OMBackend.Runtime.REDUCED_SYSTEM)
  local variableIdxToEquationIdx = G.badjlist
  local equationIdxToVariableIdx = G.fadjlist
  local idx = 0
  #= Shallow search. See if we can find the right equation directly =#
  for idx in equationIndices
    #= This equation does only depend on one variable. We are done =#
    if length(equationIdxToVariableIdx[idx]) == 1
      #= Then it is solved in this particular equation =#
      return idx
    end
    push!(usedEqIndices, idx)
  end
  for eqIdx in equationIndices
    for vIdx in equationIdxToVariableIdx[eqIdx]
      newEqIndices = filter((x) -> !(x in usedEqIndices), variableIdxToEquationIdx[vIdx])
      if isempty(newEqIndices)
        continue
      end
      idx = getRootEquation(newEqIndices; usedEqIndices = usedEqIndices)
    end
  end
  return idx
end

"""
Gets the constant value of an equation if such exist.
Throws an error otherwise
"""
function getConstantValueOfEq(eqIdx::Int, system)::Float64
  local equations = ModelingToolkit.equations(system)
  local equation = equations[eqIdx]
  @info "eq" equation
  @info "lhs" typeof(equation.lhs)
  @assert typeof(equation.lhs) != Number || typeof(equation.rhs) != Number "One side (lhs/rhs )needs to be a constant float"
  if equation.lhs isa Number
    return equation.lhs
  end
  return equation.rhs
end


end #= module =#
