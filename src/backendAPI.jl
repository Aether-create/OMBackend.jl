
#=
* This file is part of OpenModelica.
*
* Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
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

using MetaModelica
using ExportAll
using Absyn
#= For interactive evaluation. =#
using ModelingToolkit
using SymbolicUtils
using DifferentialEquations


import .Backend.BDAE
import .Backend.BDAECreate
import .Backend.BDAEUtil
import .Backend.Causalize
import ..CodeGeneration
import OMBackend.CodeGeneration
import .SimulationCode
import ..Runtime
import Base.Meta
import SCode
import JuliaFormatter
import Plots
import REPL
import OMFrontend

const latexSymbols = REPL.REPLCompletions.latex_symbols
#= Settings =#
@enum BackendMode begin
  DAE_MODE = 1
  ODE_MODE = 2 #Currently not in operation
  MTK_MODE = 3
end

function info()
  println("OMBackend.jl")
  println("A Julia backend for the Equation Oriented Language Modelica!")
  println("Run any test module by executing runExample(<model-name>)")
  println("Available Example models include:")
  for (k,v) in EXAMPLE_MODELS
    println(k)
  end
end

"""
  If this function is called a DAG representing the equations is plotted to your working directory.
"""
function plotGraph(shouldPlot::Bool)
  global PLOT_EQUATION_GRAPH = shouldPlot
end

""" 
  MTK models
"""
const COMPILED_MODELS_MTK = Dict()

function translate(frontendDAE::Union{DAE.DAE_LIST, OMFrontend.Main.FlatModel})::Tuple{String, Expr}
  bDAE = lower(frontendDAE)
  simCode = generateSimulationCode(bDAE; mode = MTK_MODE)
  return generateMTKTargetCode(simCode)
end

"""
`function dumpInitialSystem(frontendDAE::DAE.DAE_LIST)`
 Dumps a textual representation of the initial system.
"""
function dumpInitialSystem(frontendDAE::DAE.DAE_LIST)::String
  str =  "Length of frontend DAE:" * string(length(frontendDAE.elementLst)) * "\n"
  bDAE = BDAECreate.lower(frontendDAE)
  str *= BDAEUtil.stringHeading1(bDAE, "translated")
  return str
end

"""
`function printInitialSystem(frontendDAE::DAE.DAE_LIST)`
 Dumps a textual representation of the initial system.
"""
function printInitialSystem(frontendDAE::DAE.DAE_LIST)
  print(dumpInitialSystem(frontendDAE::DAE.DAE_LIST))
end

"""
 Transforms given DAE-IR/Hybrid DAE to backend DAE-IR (BDAE-IR)
"""
function lower(frontendDAE::DAE.DAE_LIST)::BDAE.BACKEND_DAE
  local bDAE::BDAE.BACKEND_DAE
  local simCode::SIM_CODE
  @debug "Length of frontend DAE:" length(frontendDAE.elementLst)
  @assert typeof(listHead(frontendDAE.elementLst)) == DAE.COMP
  #= Create Backend structure from Frontend structure =#
  bDAE = BDAECreate.lower(frontendDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "translated"));
  #= Expand arrays =#
  (bDAE, expandedVars) = Causalize.expandArrayVariables(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "Array variables expanded"));
  #= Expand Array variables in equation system=#
  bDAE = Causalize.detectAndReplaceArrayVariables(bDAE, expandedVars)
  @debug(BDAEUtil.stringHeading1(bDAE, "Equation system variables expanded"));
  #= Transform if expressions to if equations =#
  @debug(BDAEUtil.stringHeading1(bDAE, "if equations transformed"));
  bDAE = Causalize.detectIfExpressions(bDAE)
  #= Mark state variables =#
  bDAE = Causalize.detectStates(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "states marked"));
  #= We always residualize since residuals are easier to work with =#
  bDAE = Causalize.residualizeEveryEquation(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "residuals"));
  return bDAE
end


"""
  Transforms given FlatModelica to backend DAE-IR (BDAE-IR).
"""
function lower(frontendDAE::OMFrontend.Main.FlatModel)
  local bDAE = BDAECreate.lower(frontendDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "translated"));
  #= Expand arrays =#
  (bDAE, expandedVars) = Causalize.expandArrayVariables(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "Array variables expanded"));
  #= Expand Array variables in equation system=#
  bDAE = Causalize.detectAndReplaceArrayVariables(bDAE, expandedVars)
  @debug(BDAEUtil.stringHeading1(bDAE, "Equation system variables expanded"));
  #= Transform if expressions to if equations =#
  @debug(BDAEUtil.stringHeading1(bDAE, "if equations transformed"));
  bDAE = Causalize.detectIfExpressions(bDAE)
  #= Mark state variables =#
  bDAE = Causalize.detectStates(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "states marked"));
  bDAE = Causalize.residualizeEveryEquation(bDAE)
  #= =#
  @debug(BDAEUtil.stringHeading1(bDAE, "residuals"));
  return bDAE
end

"""
  Transforms  BDAE-IR to simulation code for DAE-mode
"""
function generateSimulationCode(bDAE::BDAE.BACKEND_DAE; mode)::SimulationCode.SimCode
  simCode = SimulationCode.transformToSimCode(bDAE; mode = mode)
  @debug BDAEUtil.stringHeading1(simCode, "SIM_CODE: transformed simcode")
  return simCode
end


"""
`generateMTKTargetCode(simCode::SimulationCode.SIM_CODE)`
  Generates code interfacing ModelingToolkit.jl
  The resulting code is saved in a table which contains functions that where simulated
  this session. Returns the generated modelName and corresponding generated code
"""
function generateMTKTargetCode(simCode::SimulationCode.SIM_CODE)
  #= Target code =#
  (modelName::String, modelCode::Expr) = CodeGeneration.generateMTKCode(simCode)
  @debug "Functions:" modelCode
  @debug "Model:" modelName
  modelName = replace(modelName, "." => "__")
  if haskey(COMPILED_MODELS_MTK, modelName)
    #= Check if the newly generated model is different from the previous model =#
    local changeDetected = COMPILED_MODELS_MTK[modelName] == modelCode
    COMPILED_MODELS_MTK[modelName] = (modelCode, changeDetected)
  else
    COMPILED_MODELS_MTK[modelName] = (modelCode, false)
  end
  return (modelName, modelCode)
end

function getCompiledModel(modelName)
  return first(COMPILED_MODELS_MTK[modelName])
end

"""
  Returns true if the model was compiled again
"""
function modelWasCompiledAgain(modelName)
  return last(COMPILED_MODELS_MTK[modelName])
end

"""
  Writes a model to file by default the file is formatted and comments are kept.
"""
function writeModelToFile(modelName::String, filePath::String; keepComments = true, formatFile = true)
  model = getCompiledModel(modelName)
  fileName = "$modelName.jl"
  try
    if keepComments == false
      strippedModel = CodeGeneration.stripComments(model)
    end
    local strippedModel = CodeGeneration.stripBeginBlocks(model)
    local modelStr::String = "$strippedModel"
    formattedModel = if formatFile
      JuliaFormatter.format_text(modelStr,
                                 remove_extra_newlines = true,
                                 always_use_return = true)
    else
      modelStr
    end
    writeStringToFile(filePath, formattedModel)
  catch e
    @info "Failed writing $model to file: $fileName"
    @error e
  end
end

"""
    Write the contents of a string to file.
"""
function writeStringToFile(fileName::String, contents::String)
  local fdesc = open(fileName, "w")
  write(fdesc, contents)
  close(fdesc)
end

"""
  Prints a model. 
  If the specified model exists. Print it to stdout.
"""
function printModel(modelName::String; keepComments = true, keepBeginBlocks = true)
  try
    local model::Expr
    model = getCompiledModel(modelName)
    strippedModel = "$model"
    #= Remove all the redudant blocks from the model =#
    if keepComments == false
      strippedModel = CodeGeneration.stripComments(model)
    end
    if keepBeginBlocks == false
      strippedModel = CodeGeneration.stripBeginBlocks(model)
    end    
    local modelStr::String = "$strippedModel"
    formattedResults = JuliaFormatter.format_text(modelStr;
                                                  remove_extra_newlines = true,
                                                  always_use_return = false)
    println(formattedResults)
  catch e 
    @error "Model: $(modelName) is not compiled. Available models are: $(availableModels())"
    throw("Error printing model")
  end
end

"""
    Prints available compiled models to stdout
"""
function availableModels()::String
  str = "Compiled models (MTK-MODE):\n"
    for m in keys(COMPILED_MODELS_MTK)
      str *= "  $m\n"
    end  
  return str
end

"""
`simulateModel(modelName::String; MODE = DAE_MODE ,tspan=(0.0, 1.0), solver = :solver)
  Simulates model interactivly.
The solver need to be passed with a : before the name, example:
OMBackend.simulateModel(modelName, tspan = (0.0, 1.0), solver = :(Tsit5()));
"""
function simulateModel(modelName::String; tspan=(0.0, 1.0), solver = :(Rodas5()))
  #= Strings containing . need to be in a format suitable for Julia =#
  modelName = replace(modelName, "." => "__")
  local modelCode::Expr  

  try
    modelCode = getCompiledModel(modelName)
  catch err
    println("Failed to simulate model.")
    println("Available models are:")
    availableModels()
  end
  try
    @eval $(:(import OMBackend))
    #= Below is needed to pass the custom solver=#
    strippedModel = CodeGeneration.stripBeginBlocks(modelCode)
    @eval $strippedModel
    local modelRunnable = Meta.parse("OMBackend.$(modelName)Simulate($(tspan); solver = $solver)")
    #= Run the model with the supplied tspan. =#
    @eval $modelRunnable
    #=
    The model is now compiled and a part of the OMBackend module.
    In the following path OMBackend.<modelName>Simulate
    =#
  catch err
    @info "Interactive evaluation failed: $err with mode: $(MODE)"
    @info err #TODO readd.
#      println(modelCode)
    throw(err)
  end
end

"""
  Resimulates an already compiled model given a model that is already active in th environment
  along with a set of parameters as key value pairs.
"""
function resimulateModel(modelName::String; solver = :(Rodas5()), tspan=(0.0, 1.0), parameters::Dict = Dict())
  #=
  Check if a compiled instance of the model already exists in the backend.
  If that is the case we do not have to recompile it.
  =#
  try
    local modelRunnable = Meta.parse("OMBackend.$(modelName)Simulate($(tspan); solver = $(solver))")
    @info modelRunnable
    res = @eval $modelRunnable
    return res
  catch #= The model is not compiled yet. Proceeding.. =#
    availModels = availableModels()
    @error "The model $(modelName) is not compiled.\n Available models are: $(availModels)"        
  end
end

"""
`turnOnLogging(mod = "OMBackend"::String)`\n
Turns on logging. An optional parameter `mod` can be used to specify which model should be logged for more granuality.
"""
function turnOnLogging(mod = "OMBackend"::String)
  ENV["JULIA_DEBUG"] = mod
end
