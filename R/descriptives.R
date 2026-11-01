#
# Copyright (C) 2013-2018 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

DescriptivesInternal <- function(jaspResults, dataset, options) {

  variables     <- unlist(options[["variables"]])
  variableTypes <- options[["variables.types"]]

  scaleVariables    <- variables[variableTypes == "scale"]
  nonScaleVariables <- variables[variableTypes != "scale"]

  splitName <- options$splitBy
  makeSplit <- splitName != ""
  numberMissingSplitBy <- 0

  if (makeSplit && length(variables) > 0) {
    splitFactor <- dataset[[splitName]]
    splitLevels <- levels(splitFactor)
    # remove missing values from the grouping variable
    dataset <- dataset[!is.na(splitFactor), ]

    numberMissingSplitBy <- sum(is.na(splitFactor))

    # Actually remove missing values from the split factor
    splitFactor <- na.omit(splitFactor)
    # create a list of datasets, one for each level
    # special case pareto plot could have a count variable
    if (options[["paretoAddCountVariable"]] != "") variables <- c(variables, options[["paretoAddCountVariable"]])
    splitDat <- split(dataset[variables], splitFactor)

  }


  .descriptivesDescriptivesTable(dataset, options, jaspResults, numberMissingSplitBy = numberMissingSplitBy)

  # Covariance matrix
  if (options[["covariance"]] || options[["correlation"]]) {
    if (is.null(jaspResults[["associationMatrix"]])) {
      cName <- ifelse(options[["splitBy"]] != "", "Matrices", "Matrix")
      jaspResults[["associationMatrix"]] <- createJaspContainer(gettextf("Association %s", cName))
      jaspResults[["associationMatrix"]]$dependOn(c("splitBy", "covariance", "correlation", "variables", "associationMatrixUse"))
      # jaspResults[["stemAndLeaf"]]$position <- 11
    }

    if (sum(variableTypes == "scale") > 1) {
      .descriptivesCovarianceTables(
        container = jaspResults[["associationMatrix"]],
        dataset   = if (makeSplit) splitDat else dataset,
        variables = options[["variables"]][options[["variables.types"]] == "scale"],
        options   = options
      )
    }
  }

  # Frequency table
  if (options$frequencyTables) {
    if (is.null(jaspResults[["tables"]])) {
      jaspResults[["tables"]] <- createJaspContainer(gettext("Frequency Tables"))
      jaspResults[["tables"]]$dependOn(c("frequencyTables", "splitBy", "frequencyTablesMaximumDistinctValues"))
      jaspResults[["tables"]]$position <- 3
    }

    .descriptivesFrequencyTables(dataset, options, jaspResults[["tables"]])
  }

  # Correlation plot
  if (options$correlationPlots && length(variables) > 1) {
    if (is.null(jaspResults[["matrixPlot"]])) {
      if (makeSplit) {
        jaspResults[["matrixPlot"]] <- createJaspContainer(title = gettext("Correlation plots"))
        corrPlot <- jaspResults[["matrixPlot"]]
        corrPlot$dependOn(c("correlationPlots", "splitBy", "variables"))

        for (i in seq_along(splitLevels))
          corrPlot[[splitLevels[i]]] <- .descriptivesMatrixPlot(splitDat[[i]], options, splitLevels[i])

      } else {
        jaspResults[["matrixPlot"]] <- .descriptivesMatrixPlot(dataset, options, gettext("Correlation plot")) # Create one plot
      }

      jaspResults[["matrixPlot"]]$position <- 6
    }
  }

  # Distribution plots
  if (options$distributionPlots) {
    if (is.null(jaspResults[["distributionPlots"]])) {
      jaspResults[["distributionPlots"]] <- createJaspContainer(gettext("Distribution Plots"))
      jaspResults[["distributionPlots"]]$dependOn(c(
        "variables", "distributionPlots", "splitBy", "distributionAndCorrelationPlotHistogramBinWidthType", "distributionAndCorrelationPlotDensity",
        "distributionAndCorrelationPlotRugMarks", "distributionAndCorrelationPlotHistogramManualNumberOfBins"
      ))
      jaspResults[["distributionPlots"]]$position <- 5
    }

    distPlots <- jaspResults[["distributionPlots"]]

    for (varIndex in seq_along(variables)) {
      var <- variables[varIndex]
      varType <- variableTypes[varIndex]
      if (is.null(distPlots[[var]])) {
        if (makeSplit) {
          distPlots[[var]] <- .descriptivesFrequencyPlots(dataset = splitDat, options = options,
                                                          variable = var, variableType = varType)
        } else {
          distPlots[[var]] <- .descriptivesFrequencyPlots(dataset = dataset, options = options,
                                                          variable = var, variableType = varType)
        }
      }
    }
  }

  # Box plots
  if (options$boxPlot) {
    if (is.null(jaspResults[["boxPlot"]])) {
      jaspResults[["boxPlot"]] <- createJaspContainer(gettext("Boxplots"))
      jaspResults[["boxPlot"]]$dependOn(c("boxPlot", "splitBy", "quantilesType"))
      jaspResults[["boxPlot"]]$position <- 7
    }

    splitPlots <- jaspResults[["boxPlot"]]

    for (var in scaleVariables) {
      if (is.null(splitPlots[[var]])) {
        splitPlots[[var]] <- .descriptivesSplitPlot(dataset = dataset, options = options, variable = var)
        splitPlots[[var]]$dependOn(optionContainsValue = list(variables = var))
      }
    }
  }

  # QQ plots
  if (options$qqPlot) {
    if (is.null(jaspResults[["QQPlots"]])) {
      jaspResults[["QQPlots"]] <- createJaspContainer(
        if (length(variables) > 1 || length(levels(dataset[[splitName]])) > 1) # there will be more than one Q-Q Plot
          gettext("Q-Q Plots")
        else # only one Q-Q Plot
          gettext("Q-Q Plot")
      )
      jaspResults[["QQPlots"]]$dependOn(c("qqPlot", "qqPlotCi", "qqPlotCiLevel", "splitBy", "quantilesType"))
      jaspResults[["QQPlots"]]$position <- 8
    }
    QQPlots <- jaspResults[["QQPlots"]]
    if (makeSplit) {
      qqSplitFactor <- dataset[[splitName]]
      if (length(qqSplitFactor) == 0L)
        return(createJaspPlot(error = gettext("Plotting is not possible: Variable only contains NA!"), dependencies = "splitBy"))

      # gives the different split values
      qqSplitLevels <- levels(qqSplitFactor)
      # remove missing values from the grouping variable
      dataset <- dataset[!is.na(qqSplitFactor), ]
      for (var in scaleVariables) {
        if (is.null(QQPlots[[var]])) {
            deeperQQPlots <- createJaspContainer(paste0(var))
            deeperQQPlots$dependOn(optionContainsValue = list(variables = var))
            QQPlots[[var]] <- deeperQQPlots
            # splits dataset according to split values
            qqSplitData <- split(dataset, qqSplitFactor)
            for (lev in seq_along(qqSplitLevels)) {
              QQPlots[[var]][[paste0(var, lev)]] <- .descriptivesQQPlot(dataset = qqSplitData[[lev]], options = options, qqvar = var, levelName = qqSplitLevels[lev])
            }
        }
      }
    } else { # no split
      for (var in scaleVariables) {
        if (is.null(QQPlots[[var]])) {
          QQPlots[[var]] <- .descriptivesQQPlot(dataset = dataset, options = options, qqvar = var)
        }
      }
    }
  }

  # Pie charts
  if (options$pieChart) {
    if (is.null(jaspResults[["pieCharts"]])) {
      jaspResults[["pieCharts"]] <- createJaspContainer(gettext("Pie charts"))
      jaspResults[["pieCharts"]]$dependOn(c("splitBy", "pieChart", "colorPalette"))
      jaspResults[["pieCharts"]]$position <- 9
    }

    piePlots <- jaspResults[["pieCharts"]]
    jaspGraphs::setGraphOption("palette", options[["colorPalette"]])
    for (var in nonScaleVariables) {
      if (is.null(piePlots[[var]])) {
        piePlots[[var]] <- if (makeSplit) {
          .descriptivesPieChart(dataset = splitDat, options = options, variable = var)
        } else {
          .descriptivesPieChart(dataset = dataset, options = options, variable = var)
        }
      }
    }
  }

  # Stem and leaf
  if (options[["stemAndLeaf"]]) {
    if (is.null(jaspResults[["stemAndLeaf"]])) {
      jaspResults[["stemAndLeaf"]] <- createJaspContainer(gettext("Stem and Leaf"))
      jaspResults[["stemAndLeaf"]]$dependOn(c("splitBy", "stemAndLeaf", "stemAndLeafScale"))
      jaspResults[["stemAndLeaf"]]$position <- 11
    }

    if (any(variableTypes == "scale")) {
      .descriptivesStemAndLeafTables(
        container = jaspResults[["stemAndLeaf"]],
        dataset   = if (makeSplit) splitDat else dataset,
        variables = scaleVariables,
        options   = options
      )
    }
  }

  # Scatter plots
  if (options[["scatterPlot"]]) {
    if (is.null(jaspResults[["scatterPlots"]])) {
      jaspResults[["scatterPlots"]] <- createJaspContainer(gettext("Scatter Plots"))
      jaspResults[["scatterPlots"]]$dependOn(c(
        "splitBy", "scatterPlot", "scatterPlotGraphTypeAbove", "scatterPlotGraphTypeRight", "scatterPlotRegressionLine",
        "scatterPlotRegressionLineCi", "scatterPlotRegressionLineCiLevel", "scatterPlotRegressionLineType", "scatterPlotLegend",
        "colorPalette"
      ))
      jaspResults[["scatterPlots"]]$position <- 10
    }
    .descriptivesScatterPlots(jaspResults[["scatterPlots"]], dataset, scaleVariables, splitName, options)
  }

  # Interval plots
  if (options$intervalPlot) {
    if (is.null(jaspResults[["IntervalPlots"]])) {
      jaspResults[["IntervalPlots"]] <- createJaspContainer(gettext("Interval plots"))
      jaspResults[["IntervalPlots"]]$dependOn(c("intervalPlot", "splitBy"))
      jaspResults[["IntervalPlots"]]$position <- 12
    }

    intervalPlots <- jaspResults[["IntervalPlots"]]

    for (var in scaleVariables) {
      if (is.null(intervalPlots[[var]])) {
        intervalPlots[[var]] <- .descriptivesIntervalPlot(dataset = dataset, options = options, variable = var)
        intervalPlots[[var]]$dependOn(optionContainsValue = list(variables = var))
      }
    }
  }

  # Dot plots
  if (options[["dotPlot"]]) {
    if (is.null(jaspResults[["DotPlots"]])) {
      jaspResults[["DotPlots"]] <- createJaspContainer(gettext("Dot Plots"))
      jaspResults[["DotPlots"]]$dependOn(c("splitBy", "dotPlot"))
      jaspResults[["DotPlots"]]$position <- 13
    }

    dotPlots <- jaspResults[["DotPlots"]]

    for (varIndex in seq_along(variables)) {
      var <- variables[varIndex]
      varType <- variableTypes[varIndex]
      if (is.null(dotPlots[[var]])) {
        dotPlots[[var]] <- .descriptivesDotPlots(
          dataset = if (makeSplit) splitDat else dataset,
          options = options,
          variable = var,
          variableType = varType
        )
      }
    }
  }

  # Heatmap
  if (options[["heatmapHorizontalAxis"]] != "" && options[["heatmapVerticalAxis"]] != "") {
    if (is.null(jaspResults[["heatmaps"]])) {
      jaspResults[["heatmaps"]] <- createJaspContainer(gettext("Heatmaps"))
      jaspResults[["heatmaps"]]$dependOn(c(
        "heatmapHorizontalAxis", "heatmapVerticalAxis",
        "heatmapDisplayValue", "heatmapTileWidthHeightRatio", "heatmapLegend",
        "heatmapStatisticContinuous", "heatmapStatisticDiscrete",
        "colorPalette", "splitBy", "variables", "heatmapPlot"
      ))
      jaspResults[["heatmaps"]]$position <- 14
    }

    .descriptivesHeatmaps(jaspResults[["heatmaps"]], dataset, variables, options)
  }

  # Pareto plots
  if (options[["paretoPlot"]]) {

    if (is.null(jaspResults[["paretoPlots"]])) {
      jaspResults[["paretoPlots"]] <- createJaspContainer(gettext("Pareto Plots"))
      jaspResults[["paretoPlots"]]$dependOn(c("paretoPlot", "splitBy", "paretoPlotRule", "paretoPlotRuleCi",
                                              "paretoShiftAccumulationLine", "paretoAddCountVariable",
                                              "paretoPlotTiltXAxisLabels"))
      jaspResults[["paretoPlots"]]$position <- 15
    }

    parPlots <- jaspResults[["paretoPlots"]]

    for (var in nonScaleVariables) {
      # Only categorical variables; skip the count variable — it should not get its own Pareto plot
      if (var == options[["paretoAddCountVariable"]])
        next
      if (is.null(parPlots[[var]])) {
        parPlots[[var]] <- if (makeSplit) {
          .descriptivesParetoPlots(splitDat, var, options)
        } else {
          .descriptivesParetoPlots(dataset, var, options)
        }
      }
    }
  }

  # Density plots
  if (options[["densityPlot"]]) {
    if (is.null(jaspResults[["densityPlot"]])) {
      jaspResults[["densityPlot"]] <- createJaspContainer(gettext("Frequency Plots"))
      jaspResults[["densityPlot"]]$dependOn(c(
        "densityPlot", "densityPlotSeparate", "densityPlotType", "customHistogramPosition",
        "colorPalette", "splitBy", "variables", "densityPlotTransparency"
      ))
      jaspResults[["densityPlot"]]$position <- 17
    }

    .descriptivesDensityPlots(jaspResults[["densityPlot"]], dataset, variables, options)
  }

  # Variable summary plots (nominal)
  if (options[["vSummaryPlot"]] && options[["vNominal"]]) {
    if (is.null(jaspResults[["summaryPlotNominal"]])) {
      jaspResults[["summaryPlotNominal"]] <- createJaspContainer(gettext("Variable Summary Plots (Nominal)"))
      jaspResults[["summaryPlotNominal"]]$dependOn(c(
        "vSummaryPlot", "vNominal",
        "summaryNominalType", "summaryNominalDisp", "colorPalette"
      ))
      jaspResults[["summaryPlotNominal"]]$position <- 18 # 1-17 are taken
    }

    .descriptivesSummaryNominal(jaspResults[["summaryPlotNominal"]], dataset, options)
  }

  # Variable summary plots (ordinal)
  if (options[["vSummaryPlot"]] && options[["vOrdinal"]]) {
    if (is.null(jaspResults[["summaryPlotOrdinal"]])) {
      jaspResults[["summaryPlotOrdinal"]] <- createJaspContainer(gettext("Variable Summary Plots (Ordinal)"))
      jaspResults[["summaryPlotOrdinal"]]$dependOn(c(
        "vSummaryPlot", "vOrdinal",
        "summaryOrdinalDisp", "colorPalette", "quantilesType"
      ))
      jaspResults[["summaryPlotOrdinal"]]$position <- 19
    }

    .descriptivesSummaryOrdinal(jaspResults[["summaryPlotOrdinal"]], dataset, options)
  }

  # Variable summary plots (scale)
  if (options[["vSummaryPlot"]] && options[["vScale"]]) {
    if (is.null(jaspResults[["summaryPlotScale"]])) {
      jaspResults[["summaryPlotScale"]] <- createJaspContainer(gettext("Variable Summary Plots (Scale)"))
      # the variable box is not listed: each plot depends on its own variable via optionContainsValue
      jaspResults[["summaryPlotScale"]]$dependOn(c(
        "vSummaryPlot", "vScale", "summaryScaleMainType", "summaryScaleRaindrops", "summaryScaleDensline",
        "summaryScaleQQ", "summaryScaleTimeSeries", "summaryScaleTimeSeriesDots", "summaryScaleNumbers",
        "colorPalette", "quantilesType"
      ))
      jaspResults[["summaryPlotScale"]]$position <- 20
    }

    .descriptivesSummaryScale(jaspResults[["summaryPlotScale"]], dataset, options)
  }

  # Likert plots
  if (options[["likertPlot"]]) {
    if (is.null(jaspResults[["likertPlot"]])) {
      jaspResults[["likertPlot"]] <- createJaspContainer(gettext("Likert Plots"))
      jaspResults[["likertPlot"]]$dependOn(c(
        "likertPlot", "splitBy", "variables",
        "likertPlotAssumeVariablesSameLevel", "likertPlotAdjustableFontSize"
      ))
      jaspResults[["likertPlot"]]$position <- 16
    }

    likPlots <- jaspResults[["likertPlot"]]

    for (var in nonScaleVariables) {
      # exclude non-categorical variables from dataframe
      if (is.numeric(dataset[[var]])) {
        if (makeSplit) {
          for (i in seq_along(splitLevels))
            splitDat[[i]] <- splitDat[[i]][, !names(splitDat[[i]]) %in% c(var), drop = FALSE]
        } else {
          dataset <- dataset[, !names(dataset) %in% c(var), drop = FALSE]
        }
      }
    }

    if (makeSplit) {
      for (i in seq_along(splitLevels))
        likPlots[[splitLevels[i]]] <- .descriptivesLikertPlots(splitDat[[i]], splitLevels[i], options)
    } else {
      jaspResults[["likertPlot"]] <- .descriptivesLikertPlots(dataset, gettext("Likert Plots"), options)
    }
  }

  return()
}

.descriptivesDescriptivesTable <- function(dataset, options, jaspResults, numberMissingSplitBy = 0) {
  if (!is.null(jaspResults[["stats"]])) {
    return()
  } # The options for this table didn't change so we don't need to rebuild it

  wantsSplit <- options$splitBy != ""
  variables <- unlist(options$variables)
  equalGroupsNo <- options$quantilesForEqualGroupsNumber
  percentilesPercentiles <- unique(options$percentileValues)
  stats <- createJaspTable(gettext("Descriptive Statistics"))
  stats$transpose <- !options[["descriptivesTableTransposed"]] # the table is transposed by default
  stats$position <- 1

  if (numberMissingSplitBy)
    stats$addFootnote(message = gettextf("Excluded %1$i rows from the analysis that correspond to the missing values of the split-by variable %2$s", numberMissingSplitBy, options$splitBy))

  stats$dependOn(c(
    "splitBy", "variables", "quantilesForEqualGroupsNumber", "percentileValues", "mode", "median",
    "meanArithmetic", "meanGeometric", "meanHarmonic",
    "seMean", "sd", "coefficientOfVariation", "variance", "skewness", "kurtosis", "shapiroWilkTest",
    "range", "iqr", "mad", "madRobust", "minimum", "maximum", "sum", "quartiles", "quantilesForEqualGroups",
    "percentiles", "descriptivesTableTransposed", "valid", "missing", "meanCi", "meanCiLevel", "meanCiMethod",
    "sdCi", "sdCiLevel", "sdCiMethod", "varianceCiMethod", "varianceCi", "varianceCiLevel", "ciBootstrapSamples",
    "quantilesType"
  ))

  if (wantsSplit) {
    stats$transposeWithOvertitle <- TRUE
    stats$addColumnInfo(name = "Variable", title = "", type = "string")
    stats$addColumnInfo(name = "Level",    title = "", type = "string")
  } else {
    stats$addColumnInfo(name = "Variable", title = "", type = "string")
  }

  formattedMeanCiPercent <- format(100 * options[["meanCiLevel"]], digits = 3, drop0trailing = TRUE)
  formattedSdCiPercent <- format(100 * options[["sdCiLevel"]], digits = 3, drop0trailing = TRUE)
  formattedVarianceCiPercent <- format(100 * options[["varianceCiLevel"]], digits = 3, drop0trailing = TRUE)

  # only add overtitle for CIs if table is transposed, else describe CIs in title
  if (options[["descriptivesTableTransposed"]]) {
    meanCiOvertitle <- gettextf("%s%% Confidence Interval Mean", formattedMeanCiPercent)
    sdCiOvertitle <- gettextf("%s%% Confidence Interval Std. Dev.", formattedSdCiPercent)
    varianceCiOvertitle <- gettextf("%s%% Confidence Interval Variance", formattedVarianceCiPercent)
    meanCiUbTitle <- gettext("Upper")
    meanCiLbTitle <- gettext("Lower")
    sdCiUbTitle <- gettext("Upper")
    sdCiLbTitle <- gettext("Lower")
    varianceCiUbTitle <- gettext("Upper")
    varianceCiLbTitle <- gettext("Lower")
  } else {
    meanCiOvertitle <- NULL
    sdCiOvertitle <- NULL
    varianceCiOvertitle <- NULL
    meanCiUbTitle <- gettextf("%s%% CI A. Mean Upper", formattedMeanCiPercent)
    meanCiLbTitle <- gettextf("%s%% CI A. Mean Lower", formattedMeanCiPercent)
    sdCiUbTitle <- gettextf("%s%% CI Std. Dev. Upper", formattedSdCiPercent)
    sdCiLbTitle <- gettextf("%s%% CI Std. Dev. Lower", formattedSdCiPercent)
    varianceCiUbTitle <- gettextf("%s%% CI Variance Upper", formattedVarianceCiPercent)
    varianceCiLbTitle <- gettextf("%s%% CI Variance Lower", formattedVarianceCiPercent)
  }

  if (options$valid)                          stats$addColumnInfo(name="Valid",                       title=gettext("Valid"),                   type="integer")
  if (options$missing)                        stats$addColumnInfo(name="Missing",                     title=gettext("Missing"),                 type="integer")
  if (options$mode)                           stats$addColumnInfo(name="Mode",                        title=gettext("Mode"),                    type="mixed")
  if (options$median)                         stats$addColumnInfo(name="Median",                      title=gettext("Median"),                  type="mixed")
  if (options$meanArithmetic)                 stats$addColumnInfo(name="MeanArithmetic",              title=gettext("Mean (arithmetic)"),       type="number")
  if (options$seMean)                         stats$addColumnInfo(name="Std. Error of Mean",          title=gettext("Std. Error of A. Mean"),   type="number")
  if (options$meanCi) {                       stats$addColumnInfo(name="MeanCILB",                    title=meanCiLbTitle,                      type="number", overtitle = meanCiOvertitle)
                                              stats$addColumnInfo(name="MeanCIUB",                    title=meanCiUbTitle,                      type="number", overtitle = meanCiOvertitle)}
  if (options$meanGeometric)                  stats$addColumnInfo(name="MeanGeometric",               title=gettext("Mean (geometric)"),        type="number")
  if (options$meanHarmonic)                   stats$addColumnInfo(name="MeanHarmonic",                title=gettext("Mean (harmonic)"),         type="number")
  if (options$sd)                             stats$addColumnInfo(name="Std. Deviation",              title=gettext("Std. Deviation"),          type="number")
  if (options$sdCi) {                         stats$addColumnInfo(name="SdCILB",                      title=sdCiLbTitle,                        type="number", overtitle = sdCiOvertitle)
                                              stats$addColumnInfo(name="SdCIUB",                      title=sdCiUbTitle,                        type="number", overtitle = sdCiOvertitle)}
  if (options$coefficientOfVariation)         stats$addColumnInfo(name="Coefficient of Variation",    title=gettext("Coefficient of variation"),type="number")
  if (options$mad)                            stats$addColumnInfo(name="MAD",                         title=gettext("MAD"),                     type="number")
  if (options$madRobust)                      stats$addColumnInfo(name="MAD Robust",                  title=gettext("MAD robust"),              type="number")
  if (options$iqr)                            stats$addColumnInfo(name="IQR",                         title=gettext("IQR"),                     type="number")
  if (options$variance)                       stats$addColumnInfo(name="Variance",                    title=gettext("Variance"),                type="number")
  if (options$varianceCi) {                   stats$addColumnInfo(name="VarianceCILB",                title=varianceCiLbTitle,                  type="number", overtitle = varianceCiOvertitle)
                                              stats$addColumnInfo(name="VarianceCIUB",                title=varianceCiUbTitle,                  type="number", overtitle = varianceCiOvertitle)}
  if (options$skewness) {                     stats$addColumnInfo(name="Skewness",                    title=gettext("Skewness"),                type="number")
                                              stats$addColumnInfo(name="Std. Error of Skewness",      title=gettext("Std. Error of Skewness"),  type="number") }
  if (options$kurtosis) {                     stats$addColumnInfo(name="Kurtosis",                    title=gettext("Kurtosis"),                type="number")
                                              stats$addColumnInfo(name="Std. Error of Kurtosis",      title=gettext("Std. Error of Kurtosis"),  type="number") }
  if (options$shapiroWilkTest) {              stats$addColumnInfo(name="Shapiro-Wilk",                title=gettext("Shapiro-Wilk"),            type="number")
                                              stats$addColumnInfo(name="P-value of Shapiro-Wilk",     title=gettext("P-value of Shapiro-Wilk"), type="pvalue") }
  if (options$range)                          stats$addColumnInfo(name="Range",                       title=gettext("Range"),                   type="number")
  if (options$minimum)                        stats$addColumnInfo(name="Minimum",                     title=gettext("Minimum"),                 type="mixed")
  if (options$maximum)                        stats$addColumnInfo(name="Maximum",                     title=gettext("Maximum"),                 type="mixed")

  if (options$quartiles) {
                                    stats$addColumnInfo(name="q1", title=gettext("25th percentile"), type="mixed")
                                    stats$addColumnInfo(name="q2", title=gettext("50th percentile"), type="mixed")
                                    stats$addColumnInfo(name="q3", title=gettext("75th percentile"), type="mixed")
  }

  if (options$quantilesForEqualGroups) { # I've read that there are several ways how to estimate percentiles so it should be checked if it match the SPSS way
    tempPercentNames <- 1 / equalGroupsNo * 1:(equalGroupsNo - 1) * 100

    for (i in seq_along(tempPercentNames))
      stats$addColumnInfo(name = paste0("eg", i), title = gettextf("%gth percentile", round(tempPercentNames[i], 2)), type = "mixed")
  }

  if (options$percentiles) {
    for (i in percentilesPercentiles) {
      if (i >= 0 && i <= 100) {
        stats$addColumnInfo(name = paste0("pc", i), title = gettextf("%gth percentile", i), type = "mixed")
      } else {
        .quitAnalysis(gettext("Error in Percentiles, all values should >=0 and <=100"))
      }
    }
  }

  if (options$sum)
    stats$addColumnInfo(name = "Sum", title = gettext("Sum"), type = "number")

  jaspResults[["stats"]] <- stats

  # lets just add footnotes once instead of a gazillion times..
  shouldAddNominalTextFootnote <- FALSE
  shouldAddModeMoreThanOnceFootnote <- FALSE
  shouldAddOrdinalQuantileTypeFootnote <- FALSE

  # Find the number of levels to loop over
  if (wantsSplit) {
    split <- dataset[[options$splitBy]]
    splitLevels <- levels(split)
    nLevels <- length(levels(split))

    for (varIndex in seq_along(variables)) {
      for (l in seq_len(nLevels)) {

        variable <- variables[varIndex]
        column <- dataset[[variable]][split == splitLevels[l]]
        columnType <- options[["variables.types"]][varIndex]

        subReturn <- .descriptivesDescriptivesTable_subFunction(column, columnType, list(Variable = variable, Level = splitLevels[l]), options, shouldAddNominalTextFootnote, shouldAddModeMoreThanOnceFootnote, jaspResults)

        shouldAddNominalTextFootnote                     <- subReturn$shouldAddNominalTextFootnote
        shouldAddModeMoreThanOnceFootnote                <- subReturn$shouldAddModeMoreThanOnceFootnote

        stats$addRows(subReturn$resultsCol, rowNames = paste0(variable, l))

        if (subReturn$shouldAddIdenticalFootnote) {
          stats$addFootnote(
            message = gettext("All values are identical"),
            colNames = c("Skewness", "Kurtosis", "Shapiro-Wilk", "P-value of Shapiro-Wilk"),
            rowNames = paste0(variable, l)
          )
        }

        if (subReturn$shouldAddExplainEmptySet) {
          stats$addFootnote(
            message = gettextf("Infimum (minimum) of an empty set is %1$s, supremum (maximum) of an empty set is %2$s.", "\u221E", "-\u221E"),
            colNames = c("Minimum", "Maximum"),
            rowNames = paste0(variable, l)
          )
        }

        if (subReturn$shouldAddModeMoreThanOnceFootnote)
          stats$addFootnote(message  = gettext("More than one mode exists. For nominal and ordinal data, the first mode is reported. For continuous data, the mode with the highest density estimate is reported but multiple modes may exist. We recommend visualizing the data to check for multimodality."),
                            colNames = "Mode",
                            rowNames = variable)

        if(subReturn$shouldAddGeomHarmMeansPositiveFootnote)
          stats$addFootnote(message = gettext("Geometric and harmonic means are defined only for strictly positive variables, but the data contain some non-positive values."),
                            colNames = c("MeanGeometric", "MeanHarmonic"),
                            rowNames = paste0(variable, l))

        if (!shouldAddOrdinalQuantileTypeFootnote && subReturn$shouldAddOrdinalQuantileTypeFootnote) {
          shouldAddOrdinalQuantileTypeFootnote <- TRUE
          stats$addFootnote(message = gettext("For ordinal variables, quantile type 3 was used because the selected type is not supported."))
        }

      }
    }
  } else { # we dont want to split
    for (varIndex in seq_along(variables)) {
      variable <- variables[varIndex]
      column <- dataset[[variable]]
      columnType <- options[["variables.types"]][varIndex]
      subReturn <- .descriptivesDescriptivesTable_subFunction(column, columnType, list(Variable = variable), options, shouldAddNominalTextFootnote, shouldAddModeMoreThanOnceFootnote, jaspResults)

      shouldAddNominalTextFootnote                     <- subReturn$shouldAddNominalTextFootnote
      shouldAddModeMoreThanOnceFootnote                <- subReturn$shouldAddModeMoreThanOnceFootnote

      stats$addRows(subReturn$resultsCol, rowNames = variable)

      if (subReturn$shouldAddIdenticalFootnote) {
        stats$addFootnote(
          message = gettext("All values are identical"),
          colNames = c("Skewness", "Kurtosis", "Shapiro-Wilk", "P-value of Shapiro-Wilk"),
          rowNames = variable
        )
      }

      if (subReturn$shouldAddExplainEmptySet)
        stats$addFootnote(message  = gettextf("Infimum (minimum) of an empty set is %1$s, supremum (maximum) of an empty set is %2$s.", "\u221E", "-\u221E"),
                          colNames = c("Minimum", "Maximum"),
                          rowNames = variable)

      if(subReturn$shouldAddModeMoreThanOnceFootnote)
        stats$addFootnote(message  = gettext("More than one mode exists. For nominal and ordinal data, the first mode is reported. For continuous data, the mode with the highest density estimate is reported but multiple modes may exist. We recommend visualizing the data to check for multimodality."),
                          colNames = "Mode",
                          rowNames = variable)

      if(subReturn$shouldAddGeomHarmMeansPositiveFootnote)
        stats$addFootnote(message = gettext("Geometric and harmonic means are defined only for strictly positive variables, but the data contain some non-positive values."),
                          colNames = c("MeanGeometric", "MeanHarmonic"),
                          rowNames = variable)

      if (!shouldAddOrdinalQuantileTypeFootnote && subReturn$shouldAddOrdinalQuantileTypeFootnote) {
        shouldAddOrdinalQuantileTypeFootnote <- TRUE
        stats$addFootnote(message = gettext("For ordinal variables, quantile type 3 was used because the selected type is not supported."))
      }

    }
  }

  return(stats)
}

.descriptivesDescriptivesTable_subFunction <- function(column, columnType, resultsCol, options, shouldAddNominalTextFootnote, shouldAddModeMoreThanOnceFootnote, jaspResults) {
  equalGroupsNo          <- options$quantilesForEqualGroupsNumber
  percentilesPercentiles <- unique(options$percentileValues)

  rows       <- length(column)
  na.omitted <- na.omit(column)

  # for ordinals, we show the label (a string) rather than the value
  mixedRowType <- if (columnType == "scale") "number" else "string"
  maybeToString <- if (columnType == "scale") identity else as.character
  toMixedCol <- \(x) jaspBase::createMixedColumn(maybeToString(x), mixedRowType)


  shouldAddIdenticalFootnote <- all(na.omitted[1] == na.omitted) && (options$skewness || options$kurtosis || options$shapiroWilkTest)

  valid <- length(na.omitted)
  resultsCol[["Valid"]]                   <- if (options$valid)   valid
  resultsCol[["Missing"]]                 <- if (options$missing) rows - length(na.omitted)

  if (columnType == "scale") {
    resultsCol[["Median"]]                  <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$median,              na.omitted, .descriptivesMedian, options, toMixedCol)
    resultsCol[["MeanArithmetic"]]          <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$meanArithmetic,     na.omitted, mean)
    resultsCol[["Std. Error of Mean"]]      <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$seMean,             na.omitted, function(param) { sd(param)/sqrt(length(param))} )
    resultsCol[["MeanGeometric"]]           <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$meanGeometric,      na.omitted, .geometricMean )
    resultsCol[["MeanHarmonic"]]            <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$meanHarmonic,       na.omitted, .harmonicMean )
    resultsCol[["Std. Deviation"]]          <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$sd, na.omitted, sd)
    resultsCol[["Coefficient of Variation"]]<- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$coefficientOfVariation,      na.omitted, function(param) { sd(param) / mean(param)})
    resultsCol[["MAD"]]                     <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$mad,               na.omitted, function(param) { mad(param, constant = 1) } )
    resultsCol[["MAD Robust"]]              <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$madRobust,         na.omitted, mad)
    resultsCol[["IQR"]]                     <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$iqr,               na.omitted, .descriptivesIqr, options)
    resultsCol[["Variance"]]                <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$variance,          na.omitted, var)
    resultsCol[["Kurtosis"]]                <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$kurtosis,          na.omitted, .descriptivesKurtosis)
    resultsCol[["Std. Error of Kurtosis"]]  <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$kurtosis,          na.omitted, .descriptivesSEK)
    resultsCol[["Skewness"]]                <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$skewness,          na.omitted, .descriptivesSkewness)
    resultsCol[["Std. Error of Skewness"]]  <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$skewness,          na.omitted, .descriptivesSES)
    resultsCol[["Shapiro-Wilk"]]            <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$shapiroWilkTest,   na.omitted, function(param) { res <- try(shapiro.test(param)$statistic); if(isTryError(res)) NaN else res })
    resultsCol[["P-value of Shapiro-Wilk"]] <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$shapiroWilkTest,   na.omitted, function(param) { res <- try(shapiro.test(param)$p.value);   if(isTryError(res)) NaN else res })
    resultsCol[["Sum"]]                     <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$sum,               na.omitted, sum)
    resultsCol[["Range"]]                   <- .descriptivesDescriptivesTable_subFunction_OptionChecker(options$range,             na.omitted, function(param) { range(param)[2] - range(param)[1]})
  } else if (columnType == "ordinal") {
    ordinalType <- .ordinalQuantileType(as.integer(options[["quantilesType"]]))
    if (options$median) {
      medianEst <- quantile(na.omitted, 0.5, names = FALSE, type = ordinalType)
      resultsCol[["Median"]] <- toMixedCol(medianEst)
    }
    if (options$iqr)
      resultsCol[["IQR"]] <- stats::IQR(na.omitted, type = ordinalType)
  }

  if (columnType == "scale" || columnType == "ordinal") {
    if (options$minimum) resultsCol[["Minimum"]] <- toMixedCol(min(na.omitted))
    if (options$maximum) resultsCol[["Maximum"]] <- toMixedCol(max(na.omitted))
  } else {
    if (options$minimum) resultsCol[["Minimum"]] <- toMixedCol("")
    if (options$maximum) resultsCol[["Maximum"]] <- toMixedCol("")
  }
  # validator for meanCi, sdCi, and varianceCi
  ciOptionChecker <- function(fun, na.omitted, columnType, options, jaspResults, variableName) {
    if (columnType != "scale") { # show empty cells when things cannot be computed
      return(list(upper = "", lower = ""))
    } else if (length(na.omitted) == 0L) { # show NaN when things can be computed in principle but not for this variable
      return(list(upper = NaN, lower = NaN))
    } else {
      return(fun(na.omitted, options, jaspResults, variableName))
    }
  }

  if (options[["meanCi"]]) {
    variableName <- if (is.null(resultsCol[["Level"]])) resultsCol[["Variable"]] else paste0(resultsCol[["Variable"]], resultsCol[["Level"]])
    meanCiResults <- ciOptionChecker(.descriptivesMeanCI, na.omitted, columnType, options, jaspResults, variableName)
    resultsCol[["MeanCIUB"]] <- meanCiResults$upper
    resultsCol[["MeanCILB"]] <- meanCiResults$lower
  }
  if (options[["sdCi"]]) {
    variableName <- if (is.null(resultsCol[["Level"]])) resultsCol[["Variable"]] else paste0(resultsCol[["Variable"]], resultsCol[["Level"]])
    sdCiResults <- ciOptionChecker(.descriptivesSdCI, na.omitted, columnType, options, jaspResults, variableName)
    resultsCol[["SdCIUB"]] <- sdCiResults$upper
    resultsCol[["SdCILB"]] <- sdCiResults$lower
  }
  if (options[["varianceCi"]]) {
    variableName <- if (is.null(resultsCol[["Level"]])) resultsCol[["Variable"]] else paste0(resultsCol[["Variable"]], resultsCol[["Level"]])
    varianceCiResults <- ciOptionChecker(.descriptivesVarianceCI, na.omitted, columnType, options, jaspResults, variableName)
    resultsCol[["VarianceCIUB"]] <- varianceCiResults$upper
    resultsCol[["VarianceCILB"]] <- varianceCiResults$lower
  }

  # should explain supremum and infimum of an empty set?
  shouldAddExplainEmptySet <- (options$minimum || options$maximum) && valid == 0

  shouldAddGeomHarmMeansPositiveFootnote <- columnType == "scale" && (options[["meanGeometric"]] || options[["meanHarmonic"]]) && any(na.omitted <= 0)

  shouldAddOrdinalQuantileTypeFootnote <- columnType == "ordinal" &&
    !(as.integer(options[["quantilesType"]]) %in% c(1L, 3L)) &&
    (options$median || options$iqr || options$quartiles || options$quantilesForEqualGroups || options$percentiles)

  if (options$mode) {

    if (valid == 0 || any(is.infinite(na.omitted))) {

      mode <- NA
      modeType <- "string"

    } else if (columnType == "scale") {

      # scale data -- highest density estimate
      modeType <- "number"

      if (valid == 1) {
        mode <- unique(na.omitted)
      } else {
        temp <- .descriptivesComputeModeContinuous(na.omitted)
        mode <- temp[["xValues"]][which.max(temp[["yValues"]])]

        shouldAddModeMoreThanOnceFootnote <- temp[["numModes"]] > 1L
      }

    } else {

      # ordinal or nominal data -- most frequent value
      tb <- table(na.omitted)
      allModes <- names(tb[tb == max(tb)])
      mode <- allModes[1L] # use only the first mode

      shouldAddModeMoreThanOnceFootnote <- length(allModes) > 1L
      modeType <- "string"
    }

    resultsCol[["Mode"]] <- jaspBase::createMixedRow(value = mode, type = modeType)

  }

  if (options$quartiles) {
    if (columnType == "scale" || columnType == "ordinal") {
      quartileType <- as.integer(options[["quantilesType"]])
      if (columnType == "ordinal")
        quartileType <- .ordinalQuantileType(quartileType)

      q123 <- quantile(na.omitted, c(.25, .5, .75), names = FALSE, type = quartileType)

      resultsCol[["q1"]] <- toMixedCol(q123[1])
      resultsCol[["q2"]] <- toMixedCol(q123[2])
      resultsCol[["q3"]] <- toMixedCol(q123[3])
    } else {
      resultsCol[["q1"]] <- toMixedCol("")
      resultsCol[["q2"]] <- toMixedCol("")
      resultsCol[["q3"]] <- toMixedCol("")
    }
  } else {
    resultsCol[["q1"]] <- NULL
    resultsCol[["q2"]] <- NULL
    resultsCol[["q3"]] <- NULL
  }

  equalGroupNames <- NULL

  if (options$quantilesForEqualGroups)
    equalGroupNames <- paste0("eg", seq(equalGroupsNo - 1))

  percentileNames <- NULL

  if (options$percentiles)
    percentileNames <- paste0("pc", percentilesPercentiles)

  for (row in names(resultsCol)) {
    if (startsWith(row, "eg") && ((row %in% equalGroupNames) == FALSE)) {
      resultsCol[[row]] <- NULL
    }

    if (startsWith(row, "pc") && ((row %in% percentileNames) == FALSE)) {
      resultsCol[[row]] <- NULL
    }
  }

  if (columnType == "scale" || columnType == "ordinal") {
    quartileType <- as.integer(options[["quantilesType"]])
    if (columnType == "ordinal")
      quartileType <- .ordinalQuantileType(quartileType)

    if (options$quantilesForEqualGroups) {

      for (i in seq(equalGroupsNo - 1)) {
        quantileEst <- quantile(na.omitted, c(i / equalGroupsNo), names = FALSE, type = quartileType)
        resultsCol[[paste0("eg", i)]] <- toMixedCol(quantileEst)
      }
    }

    if (options$percentiles) {

      for (i in percentilesPercentiles) {
        percentileEst <- quantile(na.omitted, c(i / 100), names = FALSE, type = quartileType)
        resultsCol[[paste0("pc", i)]] <- toMixedCol(percentileEst)
      }
    }
  } else {
    if (options$quantilesForEqualGroups) {

      for (i in seq(equalGroupsNo - 1))
        resultsCol[[paste0("eg", i)]] <- toMixedCol("")

    }

    if (options$percentiles) {

      for (i in percentilesPercentiles)
        resultsCol[[paste0("pc", i)]] <- toMixedCol("")

    }
  }

  return(list(
    resultsCol = resultsCol,
    shouldAddNominalTextFootnote = shouldAddNominalTextFootnote,
    shouldAddModeMoreThanOnceFootnote = shouldAddModeMoreThanOnceFootnote,
    shouldAddIdenticalFootnote = shouldAddIdenticalFootnote,
    shouldAddExplainEmptySet = shouldAddExplainEmptySet,
    shouldAddGeomHarmMeansPositiveFootnote = shouldAddGeomHarmMeansPositiveFootnote,
    shouldAddOrdinalQuantileTypeFootnote = shouldAddOrdinalQuantileTypeFootnote
  ))
}


.descriptivesDescriptivesTable_subFunction_OptionChecker <- function(optionToCheck, na.omitted, function_to_use, ...) {
  if (!optionToCheck)
    return(NULL)

  return(function_to_use(na.omitted, ...))
}

.descriptivesFrequencyTables <- function(dataset, options, freqTabs) {
  splitName <- options$splitBy
  wantsSplit <- splitName != ""
  splitFactor <- dataset[[splitName]]
  splitLevels <- levels(splitFactor)

  omittedVariables <- character()
  maximumDistinctValues <- options[["frequencyTablesMaximumDistinctValues"]]

  for (variableIndex in seq_along(options$variables)) {
    variable <- options[["variables"]][variableIndex]
    variableType <- options[["variables.types"]][variableIndex]
    column <- dataset[[variable]]

    noDistinctObservations <- if (variableType != "scale") nlevels(column) else length(unique(column))

    if (noDistinctObservations > maximumDistinctValues) {
      omittedVariables <- c(omittedVariables, variable)
      next
    }

    if(!is.null(freqTabs[[variable]]))
      next

    freqTab <- .descriptivesFrequencyTableMeta(variable, wantsSplit, splitName)

    freqTabs[[variable]] <- freqTab

    rows <- list()

    if (wantsSplit) {
      for (lev in splitLevels) { # also loop over the levels
        t <- table(column[splitFactor == lev])
        total <- sum(t)
        alltotal <- length(column[splitFactor == lev])
        cFreq <- 0

        for (i in seq_along(names(t))) {
          row                         <- list()
          row[["factor"]]             <- lev
          row[["Level"]]              <- names(t)[i]
          row[["Frequency"]]          <- as.vector(t[i])
          cFreq                       <- cFreq + row[["Frequency"]]
          row[["Percent"]]            <- row[["Frequency"]]/alltotal*100
          row[["Valid Percent"]]      <- row[["Frequency"]]/total*100
          row[["Cumulative Percent"]] <- cFreq/total*100
          row[[".isNewGroup"]]        <- i==1
          rows[[length(rows) + 1]]    <- row
        }

        rows[[length(rows) + 1]] <- list(
          "factor"              = "",
          "Level"               = gettext("Missing"),
          "Frequency"           = alltotal - total,
          "Percent"             = (alltotal - total) / alltotal * 100,
          "Valid Percent"       = "",
          "Cumulative Percent"  = "",
          ".isNewGroup"         = FALSE
        )

        rows[[length(rows) + 1]] <- list(
          "factor"              = "",
          "Level"               = gettext("Total"),
          "Frequency"           = alltotal,
          "Percent"             = 100,
          "Valid Percent"       = "",
          "Cumulative Percent"  = "",
          ".isNewGroup"         = FALSE
        )
      }
    } else {
      t         <- table(column)
      total     <- sum(t)
      cFreq     <- 0
      alltotal  <- length(column)

      for (lev in names(t)) {
        row                         <- list()
        row[["Level"]]              <- lev
        row[["Frequency"]]          <- as.numeric(t[lev])
        cFreq                       <- cFreq + row[["Frequency"]]
        row[["Percent"]]            <- row[["Frequency"]]/alltotal*100
        row[["Valid Percent"]]      <- row[["Frequency"]]/total*100
        row[["Cumulative Percent"]] <- cFreq/total*100
        rows[[length(rows) + 1]]    <- row
      }

      rows[[length(rows) + 1]] <- list(
        "Level"               = gettext("Missing"),
        "Frequency"           = alltotal - total,
        "Percent"             = (alltotal - total) / alltotal * 100,
        "Valid Percent"       = "",
        "Cumulative Percent"  = ""
      )

      rows[[length(rows) + 1]] <- list(
        "Level"               = gettext("Total"),
        "Frequency"           = alltotal,
        "Percent"             = 100,
        "Valid Percent"       = "",
        "Cumulative Percent"  = ""
      )
    }

    freqTab$addRows(rows)
    freqTab$status <- "complete"
  }

  if (length(omittedVariables) > 0L) {
    variableFirstTable <- setdiff(options[["variables"]], omittedVariables)
    if (length(variableFirstTable) > 0L) {
      freqTabs[["ommissionText"]] <- createJaspHtml(text =
        sprintf(
          ngettext(
            msg1 = "%2$s has more than %1$s distinct values and is omitted.",
            msg2 = "The following variables have more than %1$s distinct values and are omitted: %2$s.",
            n = length(omittedVariables)
          ),
          maximumDistinctValues,
          paste(omittedVariables, collapse = ", ")
        ), position = 1
      )
      # It's actually only the footnote that introduces this dependency...
      freqTabs$dependOn(options = "variables")
    } else {
      # all variables have more distinct values than the maximum, add a dummy table and add a footnote
      freqTab <- .descriptivesFrequencyTableMeta(NULL, wantsSplit, splitName)
      freqTab$addFootnote(gettextf(
        "All variables have more than %1$s distinct values", maximumDistinctValues
      ))
      freqTabs[["dummy"]] <- freqTab
    }
  }
}

.descriptivesFrequencyTableMeta <- function(variable = NULL, wantsSplit, splitName) {
  if (is.null(variable)) {
    freqTab <- createJaspTable(gettext("Frequencies"))
  } else {
    freqTab <- createJaspTable(gettextf("Frequencies for %s", variable))
    freqTab$dependOn(optionContainsValue = list(variables = variable))
  }

  if (wantsSplit) freqTab$addColumnInfo(name = "factor", title = splitName, type = "string", combine = TRUE)

  freqTab$addColumnInfo(name="Level",               title=variable,                       type="string")
  freqTab$addColumnInfo(name="Frequency",           title=gettext("Frequency"),           type="integer")
  freqTab$addColumnInfo(name="Percent",             title=gettext("Percent"),             type="number", format="dp:1")
  freqTab$addColumnInfo(name="Valid Percent",       title=gettext("Valid Percent"),       type="number", format="dp:1")
  freqTab$addColumnInfo(name="Cumulative Percent",  title=gettext("Cumulative Percent"),  type="number", format="dp:1")
  return(freqTab)
}

.descriptivesMatrixPlot <- function(dataset, options, name) {

  # We only make correlation plots for scale variables
  variableTypes <- unlist(options$variables.types)
  variables <- unlist(options$variables)[variableTypes == "scale"]

  l <- length(variables)
  depends <- c("correlationPlots", "variables", "variables.types", "splitBy",
               "distributionAndCorrelationPlotHistogramBinWidthType",
               "distributionAndCorrelationPlotDensity",
               "distributionAndCorrelationPlotRugMarks",
               "distributionAndCorrelationPlotHistogramManualNumberOfBins")

  if (l == 0) #Nothing to plot
    return(NULL)

  variable.statuses <- vector("list", length(variables))

  for (i in seq_along(variables)) {
    errorMessage <- .descriptivesCheckPlotErrors(dataset, variables[i], obsAmount = "< 3")
    if (!is.null(errorMessage))
      variable.statuses[[i]]$error <- errorMessage
    else
      variable.statuses[[i]]$error <- ""
  }

  plotMat <- matrix(list(), l, l)
  axisBreaks <- vector("list", l)

  # minor adjustments to plot margin to avoid cutting off the x-axis labels
  adjMargin <- ggplot2::theme(plot.margin = ggplot2::unit(c(.25, .40, .25, .25), "cm"))

  oldFontSize <- jaspGraphs::getGraphOption("fontsize")
  jaspGraphs::setGraphOption("fontsize", .85 * oldFontSize)

  # first do the diagonal and store breaks
  for (row in seq_along(variables)) {
    if (variable.statuses[[row]]$error != "") {
      plotMat[[row, row]] <- .displayErrorDescriptives(errorMessage = variable.statuses[[row]]$error)
    } else {
      plotMat[[row, row]] <- .plotMarginal(
        column         = dataset[[variables[[row]]]],
        variableName   = NULL,
        displayDensity = options[["distributionAndCorrelationPlotDensity"]],
        rugs           = options[["distributionAndCorrelationPlotRugMarks"]],
        binWidthType   = options[["distributionAndCorrelationPlotHistogramBinWidthType"]],
        numberOfBins   = options[["distributionAndCorrelationPlotHistogramManualNumberOfBins"]],
        lwd            = .7
      )
      axisBreaks[[row]] <- jaspGraphs::getAxisBreaks(plotMat[[row, row]])
    }
  }

  # now do off-diagonal and use the same breaks
  for (row in seq_len(l - 1)) {
    for (col in seq(row + 1, l)) {
      if (variable.statuses[[row]]$error != "") {
        plotMat[[row, col]] <- .displayErrorDescriptives(errorMessage = variable.statuses[[row]]$error)
      } else if (variable.statuses[[col]]$error != "") {
        plotMat[[row, col]] <- .displayErrorDescriptives(errorMessage = variable.statuses[[col]]$error)
      } else {
        plotMat[[row, col]] <- .plotScatterDescriptives(
          xVar    = dataset[[variables[[col]]]],
          yVar    = dataset[[variables[[row]]]],
          xBreaks = axisBreaks[[col]]$x,
          yBreaks = axisBreaks[[row]]$x
        ) + adjMargin
      }
    }
  }

  jaspGraphs::setGraphOption("fontsize", oldFontSize)

  # slightly adjust the positions of the labels left and above the plots.
  labelPos <- matrix(.5, 4, 2)
  labelPos[1, 1] <- .55
  labelPos[4, 2] <- .65
  p <- jaspGraphs::ggMatrixPlot(
    plotList = plotMat, leftLabels = variables, topLabels = variables,
    scaleXYlabels = NULL, labelPos = labelPos
  )

  return(createJaspPlot(plot = p, width = 250 * l + 20, aspectRatio = 1, title = name, dependencies = depends))
}

.poly.predDescriptives <- function(fit, plot = NULL, line = FALSE, xMin, xMax, lwd) {
  # create function formula
  f <- vector("character", 0)

  for (i in seq_along(coef(fit))) {
    if (i == 1) {
      temp <- paste(coef(fit)[[i]])
      f <- paste0(f, temp)
    }

    if (i > 1) {
      temp <- paste0("(", coef(fit)[[i]], ")*", "x^", i - 1)
      f <- paste(f, temp, sep = "+")
    }
  }

  x <- seq(xMin, xMax, length.out = 100)
  predY <- eval(parse(text = f))

  if (line == FALSE) {
    return(predY)
  }

  if (line) {
    plot <- plot + ggplot2::geom_line(data = data.frame(x, predY), mapping = ggplot2::aes(x = x, y = predY), size = lwd)
    return(plot)
  }
}


.plotScatterDescriptives <- function(xVar, yVar, xBreaks = NULL, yBreaks = NULL, xName = NULL, yName = NULL) {

  d <- data.frame(x = xVar, y = yVar)
  d <- na.omit(d)

  if (is.null(xBreaks))
    xBreaks <- jaspGraphs::getPrettyAxisBreaks(d$x)

  fit <- NULL

  fit <- lm(y ~ poly(x, 1, raw = TRUE), d)
  lineObj <- .poly.predDescriptives(fit, line = FALSE, xMin = xBreaks[1], xMax = xBreaks[length(xBreaks)], lwd = lwd)
  rangeLineObj <- c(lineObj[1], lineObj[length(lineObj)])
  yLimits <- range(c(pretty(yVar)), rangeLineObj)

  if (!all(is.na(yLimits))) { # this is NA in case both x and y only contain a single unique value
    if (is.null(yBreaks) || yLimits[1L] <= yBreaks[1L] || yLimits[2L] >= yBreaks[length(yBreaks)]) {
      yBreaks <- jaspGraphs::getPrettyAxisBreaks(yLimits)
    }
  }


  p <- ggplot2::ggplot(data = d, ggplot2::aes(x = x, y = y)) +
    jaspGraphs::geom_point()

  xr <- range(xBreaks)
  dfLine <- data.frame(x = xr, y = rangeLineObj)
  p <- p +
    ggplot2::geom_line(data = dfLine, ggplot2::aes(x = x, y = y), size = .7, inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(name = xName, breaks = xBreaks, limits = range(xBreaks)) +
    ggplot2::scale_y_continuous(name = yName, breaks = yBreaks, limits = range(yBreaks))

  return(jaspGraphs::themeJasp(p))
}


### empty Plot with error message ###
.displayErrorDescriptives <- function(errorMessage = NULL) {
  df <- data.frame(
    x = 0, y = 1,
    # base R version of stringr::str_wrap that automatically places \n after about 40 characters (but does not split words)
    label = paste(strwrap(errorMessage, width = 40, prefix = "\n", initial = ""), collapse = "")
  )
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y, label = label)) +
    ggplot2::geom_text(size = .4 * jaspGraphs::getGraphOption("fontsize")) +
    jaspGraphs::getEmptyTheme()
  return(p)
}


.descriptivesFrequencyPlots <- function(dataset, options, variable, variableType) {
  if (options$splitBy != "") {
    # return a collection
    split <- names(dataset)

    plotResult <- createJaspContainer(title = variable)
    plotResult$dependOn(options = "splitBy", optionContainsValue = list(variables = variable))

    for (l in split) {
      plotResult[[l]] <- .descriptivesFrequencyPlots_SubFunc(dataset = dataset[[l]],
                                                             variable = variable,
                                                             variableType = variableType,
                                                             width = options$plotWidth,
                                                             height = options$plotHeight,
                                                             displayDensity = options$distributionAndCorrelationPlotDensity,
                                                             rugs = options$distributionAndCorrelationPlotRugMarks,
                                                             title = l,
                                                             binWidthType = options$distributionAndCorrelationPlotHistogramBinWidthType,
                                                             numberOfBins = options$distributionAndCorrelationPlotHistogramManualNumberOfBins)

      plotResult[[l]]$dependOn(optionsFromObject = plotResult)
    }

    return(plotResult)
  } else {
    column <- dataset[[variable]]
    aPlot <- .descriptivesFrequencyPlots_SubFunc(dataset = dataset,
                                                 variable = variable,
                                                 variableType = variableType,
                                                 width = options$plotWidth,
                                                 height = options$plotHeight,
                                                 displayDensity = options$distributionAndCorrelationPlotDensity,
                                                 rugs = options$distributionAndCorrelationPlotRugMarks,
                                                 title = variable,
                                                 binWidthType = options$distributionAndCorrelationPlotHistogramBinWidthType,
                                                 numberOfBins = options$distributionAndCorrelationPlotHistogramManualNumberOfBins)

    aPlot$dependOn(options = "splitBy", optionContainsValue = list(variables = variable))

    return(aPlot)
  }
}

.descriptivesFrequencyPlots_SubFunc <- function(dataset, variable, variableType, width, height, displayDensity, rugs, title, binWidthType, numberOfBins) {
  freqPlot <- createJaspPlot(title = title, width = width, height = height)

  errorMessage <- .descriptivesCheckPlotErrors(dataset, variable, obsAmount = "< 3")
  column <- dataset[[variable]]
  column <- column[!is.na(column)]

  if (!is.null(errorMessage)) {
    freqPlot$setError(gettextf("Plotting not possible: %s", errorMessage))
  } else if (length(column) > 0 && variableType != "scale") {
    freqPlot$plotObject <- .barplotJASP(column, variable)
  } else if (length(column) > 0 && variableType == "scale") {
    freqPlot$plotObject <- .plotMarginal(column, variableName = variable, displayDensity = displayDensity,
                                         rugs = rugs, binWidthType = binWidthType, numberOfBins = numberOfBins)
  }

  return(freqPlot)
}

.descriptivesSplitPlot <- function(dataset, options, variable) {
  depends <- c("boxPlotColourPalette", "boxPlotViolin", "boxPlotBoxPlot", "boxPlotJitter", "boxPlotOutlierLabel")

  thePlot <- createJaspPlot(title = variable, width = options$plotWidth, height = options$plotHeight, dependencies = depends)

  additionalChecks <- if (options[["boxPlotViolin"]]) "variance"  else NULL
  minObsAmount <-     if (options[["boxPlotViolin"]]) "< 2"       else "< 1"
  errorMessage <- .descriptivesCheckPlotErrors(dataset, variable,
                                               obsAmount = minObsAmount,
                                               obsGrouping = options[["splitBy"]],
                                               additionalChecks = additionalChecks)

  if (!is.null(errorMessage)) {
    thePlot$setError(gettextf("Plotting not possible: %s", errorMessage))
  } else if (!(options$boxPlotViolin || options$boxPlotBoxPlot || options$boxPlotJitter)) {
    thePlot$setError(gettext("Plotting is not possible: No plot type selected!"))
  } else {
    # we need to know which index in y is related to which index in the actual data, so we should not forget the NAs somehow, lets make a list of indices.
    yWithNA         <- dataset[[variable]]
    y               <- na.omit(dataset[[variable]])
    yIndexToActual  <- y
    yWithNAIndex    <- 1
    yNoNAIndex      <- 1

    while (yWithNAIndex <= length(yWithNA)) {
      if (!is.na(yWithNA[[yWithNAIndex]])) {
        yIndexToActual[[yNoNAIndex]] <- row.names(dataset)[[yWithNAIndex]]
        yNoNAIndex <- yNoNAIndex + 1
      }

      yWithNAIndex <- yWithNAIndex + 1
    }

    if (is.null(dataset[[options$splitBy]])) {
      group     <- factor(rep("", length(y)))
      xlab      <- "Total"
      boxWidth  <- 0.2
      vioWidth  <- 0.3
    } else {
      group     <- dataset[[options$splitBy]][!is.na(dataset[[variable]])]
      xlab      <- options$splitBy
      boxWidth  <- 0.4
      vioWidth  <- 0.6
    }

    plotDat <- data.frame(group = group, y = y)
    row.names(plotDat) <- yIndexToActual

    # Identify outliers to label. Note that ggplot uses the unchangeable quantiles(type=7),
    # if we ever change the quantile type then the boxplot needs to be overwritten with stat_summary(geom='boxplot')
    plotDat$outlier <- FALSE

    for (level in levels(plotDat$group)) {
      v <- plotDat[plotDat$group == level, ]$y
      quantiles <- quantile(v, probs = c(0.25, 0.75), type = as.integer(options[["quantilesType"]]))
      obsIQR <- quantiles[2] - quantiles[1]
      plotDat[plotDat$group == level, ]$outlier <- v < (quantiles[1] - 1.5 * obsIQR) | v > (quantiles[2] + 1.5 * obsIQR)
    }

    plotDat$label <- ifelse(plotDat$outlier, row.names(plotDat), "")

    if (options[["boxPlotColourPalette"]]) {
      thePlot$dependOn("colorPalette") # only add color as dependency if the user wants it
      palette <- options[["colorPalette"]]
      p <- ggplot2::ggplot(plotDat, ggplot2::aes(x = group, y, fill = group)) +
        jaspGraphs::scale_JASPfill_discrete(palette) +
        jaspGraphs::scale_JASPcolor_discrete(palette)
    } else {
      p <- ggplot2::ggplot(plotDat, ggplot2::aes(x = group, y, fill = group)) +
        ggplot2::scale_fill_manual(values = rep("grey", nlevels(group))) +
        ggplot2::scale_colour_manual(values = rep("grey", nlevels(group)))
    }

    if (options[["boxPlotViolin"]]) {
      p <- p + ggplot2::geom_violin(trim = FALSE, size = 0.75, scale = "width", width = vioWidth)
    }

    if (options[["boxPlotBoxPlot"]]) {
      # if we add jittered data points, don't show outlier dots
      outlierShape <- if (options[["boxPlotJitter"]]) NA else 19

      # create summaries per group for boxplots
      # needed to make the quantiles based on the type specified by the user
      boxData <- dplyr::summarise(
        dplyr::group_by(plotDat, group),
        ymin = min(y[!outlier]), # lower end of whiskers smallest observation that is not an outlier
        lower = quantile(y, 0.25, type = as.integer(options[["quantilesType"]])),
        middle = quantile(y, 0.50, type = as.integer(options[["quantilesType"]])),
        upper = quantile(y, 0.75, type = as.integer(options[["quantilesType"]])),
        ymax = max(y[!outlier]) # upper end of whiskers largest observation that is not an outlier
      )

      # plot boxplots with errorbars based on the computed statistics
      p <- p +
        ggplot2::geom_errorbar(data = boxData, ggplot2::aes(
          ymin = ymin,
          ymax = ymax,
          x = group
        ),
        inherit.aes = FALSE, linewidth = 0.75, width = boxWidth / 2) +

        ggplot2::geom_boxplot(data = boxData, ggplot2::aes(
          x = group,
          fill = group,
          ymin = ymin,
          lower = lower,
          middle = middle,
          upper = upper,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        size = 0.75, width = boxWidth,
        stat = "identity")

      # add outliers if there are any
      # since boxplot function cannot add them due to precomputed stats
      if(any(plotDat$outlier)) {
        p <- p +
          ggplot2::geom_point(data = plotDat[plotDat$outlier, ],
                              shape = outlierShape,
                              size = 2)
      }
    }

    if (options[["boxPlotJitter"]]) {
      p <- p + ggplot2::geom_jitter(ggplot2::aes(fill = group, alpha = 0.5), size = 3, shape = 21, stroke = 1, position = ggplot2::position_jitter(width = 0.05, height = 0))
      if (options[["boxPlotColourPalette"]])
        p <- p + jaspGraphs::scale_JASPfill_discrete(palette)
      else
        p <- p + ggplot2::scale_fill_manual(values = rep("grey", nlevels(group)))
    }

    if (options[["boxPlotOutlierLabel"]]) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = label), hjust = -1, min.segment.length = 99)
    }

    ### Theming & Cleaning
    yBreaks <- if (options[["boxPlotViolin"]]) {
      jaspGraphs::getPrettyAxisBreaks(range(unlist(tapply(y, group, function(x) range(x, density(x)$x)), use.names = FALSE)))
    } else {
      jaspGraphs::getPrettyAxisBreaks(p[["data"]][["y"]])
    }
    yLimits <- range(yBreaks)

    p <- p +
      ggplot2::xlab(xlab) +
      ggplot2::scale_y_continuous(name = variable, breaks = yBreaks, limits = yLimits) +
      jaspGraphs::geom_rangeframe(sides = "l") +
      jaspGraphs::themeJaspRaw()

    thePlot$plotObject <- p
  }
  return(thePlot)
}

.descriptivesIntervalPlot <- function(dataset, options, variable) {
  thePlot <- createJaspPlot(title = variable)

  errorMessage <- .descriptivesCheckPlotErrors(dataset, variable, obsAmount = "< 1")
  if (!is.null(errorMessage)) {
    thePlot$setError(gettextf("Plotting not possible: %s", errorMessage))
    return(thePlot)
  } else {
    y <- na.omit(dataset[[variable]])
  }

  if (is.null(dataset[[options$splitBy]])) {
    group <- factor(rep("", length(y)))
    xlab <- gettext("Total")
  } else {
    group <- dataset[[options$splitBy]][!is.na(dataset[[variable]])]
    xlab <- options$splitBy
  }

  plotDat <- data.frame(group = group, y = y)

  cdata <- aggregate(y ~ group, data = plotDat, function(x) {
    mu <- mean(x)
    err <- qnorm(0.975) * (sd(x) / sqrt(length(x)))
    c(
      mean = mu,
      lower = mu - err,
      upper = mu + err
    )
  })
  df <- cbind.data.frame(group = factor(cdata[["group"]]), cdata[["y"]])

  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(df[["lower"]], df[["upper"]]))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = mean, ymin = lower, ymax = upper)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(width = .1) +
    ggplot2::xlab(xlab) +
    ggplot2::scale_y_continuous(name = variable, breaks = yBreaks, limits = range(yBreaks)) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()

  thePlot$plotObject <- p

  return(thePlot)
}

.plotMarginal <- function(column, variableName,
                          rugs = FALSE, displayDensity = FALSE, binWidthType = c("doane", "fd", "scott", "sturges", "manual"),
                          numberOfBins = NA,
                          lwd = 1) {
  binWidthType <- match.arg(binWidthType)
  variable <- na.omit(column)

  if (length(variable) == 0)
    return(NULL)

  if (binWidthType == "doane") { # https://en.wikipedia.org/wiki/Histogram#Doane's_formula
    sigma.g1 <- sqrt((6 * (length(variable) - 2)) / ((length(variable) + 1) * (length(variable) + 3)))
    g1 <- mean(abs(variable)^3)
    k <- 1 + log2(length(variable)) + log2(1 + (g1 / sigma.g1))
    binWidthType <- k
  } else if (binWidthType == "fd" && nclass.FD(variable) > 10000) { # FD-method will produce extreme number of bins and crash ggplot, mention this in footnote
    binWidthType <- 10000
  } else if (binWidthType == "manual") {
    binWidthType <- numberOfBins
  }


  h <- hist(variable, plot = FALSE, breaks = binWidthType)

  if (!displayDensity) {
    yhigh <- max(h$counts)
  } else {
    dens <- density(variable)
    yhigh <- max(max(h$density), max(dens$y))
  }

  xticks <- base::pretty(c(variable, h$breaks), min.n = 3)

  if (!displayDensity) {
    p <-
      jaspGraphs::drawAxis(
        xName = variableName, yName = gettext("Counts"), xBreaks = xticks,
        yBreaks = base::pretty(c(0, h$counts)), force = TRUE, xLabels = xticks
      )
  } else {
    p <-
      jaspGraphs::drawAxis(
        xName = variableName, yName = gettext("Density"), xBreaks = xticks,
        yBreaks = c(0, 1.05 * yhigh), force = TRUE, yLabels = NULL,
        xLabels = xticks
      )
  }


  if (displayDensity) {
    p <- p +
      ggplot2::geom_histogram(
        data    = data.frame(variable),
        mapping = ggplot2::aes(x = variable, y = ..density..),
        breaks  = h[["breaks"]],
        fill    = "grey",
        col     = "black",
        size    = .7
      ) +
      ggplot2::geom_line(
        data    = data.frame(x = dens$x, y = dens$y),
        mapping = ggplot2::aes(x = x, y = y),
        lwd     = lwd,
        col     = "black"
      )
  } else {
    p <- p +
      ggplot2::geom_histogram(
        data     = data.frame(variable),
        mapping  = ggplot2::aes(x = variable, y = ..count..),
        breaks   = h[["breaks"]],
        fill     = "grey",
        col      = "black",
        size     = .7
      )
  }

  if (rugs)
    p <- p + ggplot2::geom_rug(data = data.frame(variable), mapping = ggplot2::aes(x = variable), sides = "b")

  # JASP theme
  p <- jaspGraphs::themeJasp(p,
    axisTickWidth = .7,
    bty = list(type = "n", ldwX = .7, lwdY = 1)
  )
  # TODO: Fix jaspgraphs axis width X vs Y. See @vandenman.

  if (displayDensity)
    p <- p + ggplot2::theme(axis.ticks.y = ggplot2::element_blank())

  return(p)
}

.descriptivesDotPlots <- function(dataset, options, variable, variableType) {
  if (options$splitBy != "") {
    # return a collection
    levelsSplitFactor <- names(dataset)

    plotContainer <- createJaspContainer(title = variable)
    plotContainer$dependOn(optionContainsValue = list(variables = variable))

    for (level in levelsSplitFactor)
      plotContainer[[level]] <- .descriptivesDotPlots_SubFunc(dataset = dataset[[level]],
                                                              variable = variable,
                                                              variableType = variableType,
                                                              title = level)

    return(plotContainer)
  } else {
    dotplot <- .descriptivesDotPlots_SubFunc(dataset = dataset,
                                             variable = variable,
                                             variableType = variableType,
                                             title = variable)
    dotplot$dependOn(optionContainsValue = list(variables = variable))

    return(dotplot)
  }
}

.descriptivesDotPlots_SubFunc <- function(dataset, variable, variableType, title) {
  dotPlot <- createJaspPlot(title = title)
  x <- na.omit(dataset[[variable]])
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    dotPlot$setError(gettext("No non-missing values!"))
    return(dotPlot)
  }

  dotsize <- 1

  if (variableType == "nominal" || variableType == "ordinal") {
    tb <- as.data.frame(table(x))
    scaleX <- ggplot2::scale_x_discrete(limits = factor(tb[, 1L]))

    # so the geom_dotplot is very weird and the dot sizes are independent of the y-axis.
    # hence to avoid exceeding the y-axis, we need to make the dots smaller...
    # see also the many issues on the ggplot repo about the dotplot...
    # for example https://github.com/tidyverse/ggplot2/issues/2203
    # this post provides alternatives we should consider for the current implementation
    # https://stackoverflow.com/questions/53697235/ggplot-dotplot-what-is-the-proper-use-of-geom-dotplot
    maxFreq <- max(tb[["Freq"]])
    if (maxFreq >= 17L)
      dotsize <- 17 / maxFreq

  } else if (variableType == "scale") {

    xBreaks <- jaspGraphs::getPrettyAxisBreaks(x)
    scaleX <- ggplot2::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks))

    if (length(unique(x)) == 1)
      dotsize <- .03
  }

  p <- ggplot2::ggplot(data = data.frame(x = x), ggplot2::aes(x = x)) +
    ggplot2::geom_dotplot(binaxis = "x", stackdir = "up", fill = "grey", dotsize = dotsize) +
    ggplot2::xlab(variable) +
    ggplot2::ylab(NULL) +
    scaleX +
    jaspGraphs::geom_rangeframe(sides = "b") +
    jaspGraphs::themeJaspRaw() +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )

  dotPlot$plotObject <- p

  return(dotPlot)
}

.barplotJASP <- function(column, variable) {
  tb <- as.data.frame(table(column))
  p <- ggplot2::ggplot(data = data.frame(x = tb[, 1], y = tb[, 2]), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat = "identity", fill = "grey", col = "black", size = .3) +
    ggplot2::xlab(variable) +
    ggplot2::ylab(gettext("Counts")) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw() +
    ggplot2::theme(plot.margin = ggplot2::margin(5))

  return(p)
}

.descriptivesKurtosis <- function(x) {
  # Kurtosis function as in SPSS:
  # http://www.ats.ucla.edu/stat/mult_pkg/faq/general/kurtosis.htm
  # http://en.wikipedia.org/wiki/Kurtosis#Estimators_of_population_kurtosis

  n <- length(x)
  s4 <- sum((x - mean(x))^4)
  s2 <- sum((x - mean(x))^2)
  v <- s2 / (n - 1)
  a <- (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))
  b <- s4 / (v^2)
  c <- (-3 * (n - 1)^2) / ((n - 2) * (n - 3))

  kurtosis <- a * b + c

  return(kurtosis)
}

# R's quantile() on ordered factors only supports types 1 and 3.
# Fall back to type 3 when the user-selected type is incompatible.
.ordinalQuantileType <- function(type) {
  if (type %in% c(1L, 3L)) type else 3L
}

.descriptivesMedian <- function(x, options, toMixedCol) {
  type <- as.integer(options[["quantilesType"]])
  return(toMixedCol(quantile(x, 0.5, names = FALSE, type = type)))
}

.descriptivesIqr <- function(x, options) {
  # Interquartile range based on the stats package using specified quantile type
  type <- as.integer(options[["quantilesType"]])
  return(stats::IQR(x, type = type))
}

.descriptivesSkewness <- function(x) {
  # Skewness function as in SPSS (for samlpes spaces):
  # http://suite101.com/article/skew-and-how-skewness-is-calculated-in-statistical-software-a231005

  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  z <- (x - m) / s # z scores
  a <- n / ((n - 1) * (n - 2))

  skewness <- sum(z^3) * a

  return(skewness)
}

.descriptivesSES <- function(x) {
  # Standard Error of Skewness
  # Formula found http://web.ipac.caltech.edu/staff/fmasci/home/statistics_refs/SkewStatSignif.pdf

  n <- length(x)
  SES <- sqrt((6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3))))

  return(SES)
}

.descriptivesSEK <- function(x) {
  # Standard Error of Kurtosis
  # Formula found http://web.ipac.caltech.edu/staff/fmasci/home/statistics_refs/SkewStatSignif.pdf

  n <- length(x)
  SEK <- 2 * .descriptivesSES(x) * sqrt((n^2 - 1) / ((n - 3) * (n + 5)))

  return(SEK)
}

.descriptivesMeanCI <- function(data, options, jaspResults, variableName) {
  ciWidth <- options[["meanCiLevel"]]
  if (options[["meanCiMethod"]] == "normalModel") {
    xBar <- mean(data)
    se <- sd(data) / sqrt(length(data))
    z <- qnorm((1 - ciWidth) / 2, lower.tail = FALSE)
    errorMargin <- abs(z * se)
    CIs <- c(xBar - errorMargin, xBar + errorMargin)
  } else if (options[["meanCiMethod"]] == "bootstrap") {
    stateContainerName <- paste0("bootstrapSamples", variableName)
    means <- .bootstrapStats(data, options, jaspResults, stateContainerName)$means
    percentiles <- (1 + c(-ciWidth, ciWidth)) / 2
    CIs <- quantile(means, probs = percentiles)
  } else if (options[["meanCiMethod"]] == "oneSampleTTest" && (length(unique(data)) > 2)) {
    ttestResult <- stats::t.test(x = data, conf.level = ciWidth)
    CIs <- ttestResult[["conf.int"]]
  } else {
    CIs <- c(NA, NA)
  }
  return(list(
    "lower" = CIs[1],
    "upper" = CIs[2]
  ))
}

# Formula for chi-square confidence intervals for std. dev. and variance found here: https://en.wikipedia.org/wiki/Standard_deviation#Confidence_interval_of_a_sampled_standard_deviation
.descriptivesSdCI <- function(data, options, jaspResults, variableName) {
  ciWidth <- options[["sdCiLevel"]]
  if (options[["sdCiMethod"]] == "bootstrap") {
    stateContainerName <- paste0("bootstrapSamples", variableName)
    sds <- .bootstrapStats(data, options, jaspResults, stateContainerName)$sds
    percentiles <- (1 + c(-ciWidth, ciWidth)) / 2
    CIs <- stats::quantile(sds, probs = percentiles)
    lowerBound <- CIs[1]
    upperBound <- CIs[2]
  } else if (options[["sdCiMethod"]] == "chiSquaredModel") {
    sigma <- stats::sd(data)
    n <- length(data)
    alpha <- 1 - ciWidth
    chiSqLower <- stats::qchisq(alpha / 2, df = n - 1)
    chiSqUpper <- stats::qchisq(1 - alpha / 2, df = n - 1)
    lowerBound <- sqrt((n - 1) * sigma^2 / chiSqUpper)
    upperBound <- sqrt((n - 1) * sigma^2 / chiSqLower)
  }
  return(list(
    "upper" = upperBound,
    "lower" = lowerBound
  ))
}

.descriptivesVarianceCI <- function(data, options, jaspResults, variableName) {
  ciWidth <- options[["varianceCiLevel"]]
  if (options[["varianceCiMethod"]] == "bootstrap") {
    stateContainerName <- paste0("bootstrapSamples", variableName)
    variances <- .bootstrapStats(data, options, jaspResults, stateContainerName)$variances
    percentiles <- (1 + c(-ciWidth, ciWidth)) / 2
    CIs <- stats::quantile(variances, probs = percentiles)
    lowerBound <- CIs[1]
    upperBound <- CIs[2]
  } else if (options[["varianceCiMethod"]] == "chiSquaredModel") {
    sigma <- stats::sd(data)
    n <- length(data)
    alpha <- 1 - ciWidth
    chiSqLower <- stats::qchisq(alpha / 2, df = n - 1)
    chiSqUpper <- stats::qchisq(1 - alpha / 2, df = n - 1)
    lowerBound <- (n - 1) * sigma^2 / chiSqUpper
    upperBound <- (n - 1) * sigma^2 / chiSqLower
  }
  return(list(
    "upper" = upperBound,
    "lower" = lowerBound
  ))
}

.bootstrapStats <- function(data, options, jaspResults, stateContainerName) {
  if (!is.null(jaspResults[[stateContainerName]]$object))
    return(jaspResults[[stateContainerName]]$object)

  bootstrapSamples <- createJaspState()
  k <- options[["ciBootstrapSamples"]]
  means <- numeric(k)
  sds <- numeric(k)
  variances <- numeric(k)
  n <- length(data)
  for (i in seq_len(k)) {
    bootData <- sample(data, size = n, replace = TRUE)
    means[i] <- mean(bootData)
    sds[i] <- sd(bootData)
    variances[i] <- var(bootData)
  }
  bootstrapSamples$object <- list(means = means, sds = sds, variances = variances)
  jaspResults[[stateContainerName]] <- bootstrapSamples
  jaspResults[[stateContainerName]]$dependOn(options = c("ciBootstrapSamples"))
  return(list(means = means, sds = sds, variances = variances))
}

.descriptivesComputeModeContinuous <- function(x) {

  n <- 2^15
  bw <- stats::bw.nrd0(x)
  lowsup <- -Inf
  uppsup <- Inf
  # based on multimode::nmodes
  fn <- stats::density(x, bw = bw, n = n)
  z  <- c(1:(n - 1))
  re <- z[diff(fn$y) > 0]
  z2 <- c(1:length(re))
  se <- z2[diff(re) > 1]
  posic <- re[se]
  if (re[length(re)] < (n - 1)) {
    posic <- c(posic, re[length(re)])
  }
  posic <- posic[fn$x[posic] > lowsup]
  posic <- posic[fn$x[posic] < uppsup]
  num  <- length(posic)
  # end of multimode::nmodes
  xValues <- fn[["x"]][posic]
  yValues <- fn[["y"]][posic]
  return(list(
    numModes = num,
    xValues  = xValues,
    yValues  = yValues,
    bw       = bw
  ))
}


.descriptivesCovarianceTables <- function(container, dataset, variables, options) {

  # only make tables for which the checkbox is checked:
  for (thisAss in  c("Covariance", "Correlation")[c(options[["covariance"]], options[["correlation"]])]) {

    container[[thisAss]] <- createJaspContainer(title = thisAss)

    if (!is.data.frame(dataset)) { # dataset is split
      splitLevels <- names(dataset)
      for (split in splitLevels) {
        tableName <- paste0(thisAss, "_", split)
        if (is.null(container[[thisAss]][[tableName]])) {
          container[[thisAss]][[tableName]] <- .descriptivesCovarianceCreateSingleTable(dataset[[split]][variables], thisAss, options[["associationMatrixUse"]], split)
        }
      }
    } else {
      container[[thisAss]] <- .descriptivesCovarianceCreateSingleTable(dataset[variables], thisAss, options[["associationMatrixUse"]], thisAss)
    }
  }
}

.descriptivesCovarianceCreateSingleTable <- function(x, thisAss, useObs, title = "") {
  tb <- createJaspTable(title = title)

  tb$addColumnInfo(name = "colnames(x)", title = "", type = "string")
  for (var in colnames(x)) {
    tb$addColumnInfo(name = var, type = "number")
  }

  result <- try({
    if (thisAss == "Covariance") {
      cov(x, use = useObs)
    } else if (thisAss == "Correlation") {
      cor(x, use = useObs)
    }
  })

  if (isTryError(result)) {
    if (grepl(pattern = "missing observations", result[[1]])) {
      tb$setError(gettext("Association cannot be computed due to missing values. Try a different 'Use' setting."))
    } else {
      tb$setError(result[[1]])
    }
  }

  if (anyNA(result)) {
    tb$addFootnote(message = gettext("Some associations cannot be computed due to missing values. Try a different 'Use' setting."))
  }

  tb$setData(cbind(colnames(x), as.data.frame(result)))

  return(tb)
}


.descriptivesQQPlot <- function(dataset, options, qqvar, levelName = NULL) {
  # to put a subtitle if there is a split
  title <- qqvar
  if(!is.null(levelName))
    title <- levelName

  descriptivesQQPlot <- createJaspPlot(width = 400, aspectRatio = 1, title = title)

  if (!is.null(qqvar)) {
    errorMessage <- .descriptivesCheckPlotErrors(dataset, qqvar, obsAmount = "< 1")
    if (!is.null(errorMessage)) {
      descriptivesQQPlot$setError(gettextf("Plotting not possible: %s", errorMessage))
    } else if (sum(!is.na(dataset[[qqvar]])) < 2) {
      descriptivesQQPlot$setError(gettextf("Plotting not possible: Variable '%s' has fewer than 2 valid (non-missing) observations.", qqvar))
    } else {

      ciLevel <- if (options[["qqPlotCi"]])  options[["qqPlotCiLevel"]] else NULL

      descriptivesQQPlot$plotObject <- .descriptivesQQNormPlot(dataset[[qqvar]],
                                                               ciLevel = ciLevel,
                                                               yName = gettext("Standardized sample quantiles"),
                                                               identicalAxes = TRUE)
    }
  }

  if (is.null(levelName))
    descriptivesQQPlot$dependOn(optionContainsValue=list(variables=qqvar))
  return(descriptivesQQPlot)
}

.descriptivesPieChart <- function(dataset, options, variable) {
  if (options$splitBy != "") {
    # return a collection
    split <- names(dataset)

    plotResult <- createJaspContainer(title = variable)
    plotResult$dependOn(optionContainsValue = list(variables = variable))

    for (l in split) {
      plotResult[[l]] <- .descriptivesPieChart_SubFunc(
        dataset = dataset[[l]], variable = variable, width = options$plotWidth, height = options$plotHeight, title = l,
        palette = options[["colorPalette"]]
      )
      plotResult[[l]]$dependOn(optionsFromObject = plotResult)
    }

    return(plotResult)
  } else {
    aPlot <- .descriptivesPieChart_SubFunc(
      dataset = dataset, variable = variable, width = options$plotWidth, height = options$plotHeight, title = variable,
      palette = options[["colorPalette"]]
    )
    aPlot$dependOn(options = "splitBy", optionContainsValue = list(variables = variable))

    return(aPlot)
  }
}

.descriptivesPieChart_SubFunc <- function(dataset, variable, width, height, title, palette) {
  pieChart <- createJaspPlot(title = title, width = width, height = height)

  errorMessage <- .descriptivesCheckPlotErrors(dataset, variable, obsAmount = "< 3")
  column <- dataset[[variable]]
  column <- column[!is.na(column)]
  if (!is.null(errorMessage)) {
    pieChart$setError(gettextf("Plotting not possible: %s", errorMessage))
  } else if (length(column) > 0) {
    tb <- as.data.frame(table(column))
    pieChart$plotObject <- jaspGraphs::plotPieChart(tb[, 2], tb[, 1],
      legendName = variable,
      palette = palette
    )
  }

  return(pieChart)
}

.descriptivesScatterPlots <- function(jaspContainer, dataset, variables, split, options, name = NULL, dependOnVariables = TRUE) {
  jaspGraphs::setGraphOption("palette", options[["colorPalette"]])

  if (!is.null(split) && split != "") {
    # group with split will be defined later, to handle missingness
    legendTitle <- split
  } else {
    group <- NULL
    split <- NULL
    legendTitle <- NULL
  }

  nvar <- length(variables)

  if (nvar < 2L)
    return()

  for (i in 1:(nvar - 1L)) {
    for (j in (i + 1L):nvar) {
      v1 <- variables[i]
      v2 <- variables[j]

      # omit NAs, including split if relevant
      pairwiseData <- na.omit(dataset[, c(v1, v2, split)])
      if (!is.null(split) && split != "") {
        group <- pairwiseData[[split]]
      }

      if (!is.null(name)) {
        plotName <- name
      } else {
        plotName <- paste(v1, "-", v2)
      }

      if (is.null(jaspContainer[[plotName]])) {
        scatterPlot <- createJaspPlot(title = plotName)
        if (dependOnVariables) {
          scatterPlot$dependOn(optionContainsValue = list(variables = c(v1, v2)))
        }

        wantsAnyDensity <- any(options[c("scatterPlotGraphTypeAbove", "scatterPlotGraphTypeRight")] == "density")
        obsGrouping <- if (wantsAnyDensity) options[["splitBy"]] else NULL
        errorMessage <- .descriptivesCheckPlotErrors(pairwiseData, c(v1, v2, split),
                                                     obsAmount = "< 2",
                                                     obsGrouping = obsGrouping)
        if (is.null(errorMessage)) {

          p <- try(jaspGraphs::JASPScatterPlot(
            x                 = pairwiseData[[v1]],
            y                 = pairwiseData[[v2]],
            group             = group,
            xName             = v1,
            yName             = v2,
            showLegend        = options[["scatterPlotLegend"]],
            addSmooth         = options[["scatterPlotRegressionLine"]],
            addSmoothCI       = options[["scatterPlotRegressionLineCi"]],
            smoothCIValue     = options[["scatterPlotRegressionLineCiLevel"]],
            forceLinearSmooth = options[["scatterPlotRegressionLineType"]] == "linear",
            plotAbove         = options[["scatterPlotGraphTypeAbove"]],
            plotRight         = options[["scatterPlotGraphTypeRight"]],
            legendTitle       = legendTitle
          ))

          if (isTryError(p)) {
            errorMessage <- .extractErrorMessage(p)
          }
        }

        if (!is.null(errorMessage)) {
          scatterPlot$setError(gettextf("Plotting not possible: %s", errorMessage))
        } else {
          scatterPlot$plotObject <- p
        }

        jaspContainer[[plotName]] <- scatterPlot
      }
    }
  }
}

.descriptivesCheckPlotErrors <- function(dataset, vars, obsAmount, additionalChecks = NULL, obsGrouping = NULL) {
  errors <- .hasErrors(dataset, all.target = vars, message = "short", observations.grouping = obsGrouping,
                       type = c("infinity", "observations", additionalChecks),
                       observations.amount = obsAmount)
  if (!isFALSE(errors)) {
    return(errors$message)
  }

  return(NULL)
}

.descriptivesStemAndLeafTables <- function(container, dataset, variables, options, width = 120, atom = 1e-08) {
  # parameters for graphics::stem
  scale <- options[["stemAndLeafScale"]]
  # width and atom are set to R's default values. We could allow the user to set them but the R documentation is unclear as to what it does.

  if (!is.data.frame(dataset)) { # dataset is split

    splitLevels <- names(dataset)

    for (var in variables) {
      if (is.null(container[[var]])) {
        subcontainer <- createJaspContainer(title = var)
        subcontainer$dependOn(optionContainsValue = list(variables = var))
        container[[var]] <- subcontainer

        for (split in splitLevels) {
          tableName <- paste0("stem_and_leaf_", var, "_", split)
          subcontainer[[tableName]] <- .descriptivesStemAndLeafCreateSingleTable(dataset[[split]][[var]], split, scale, width, atom)
        }
      }
    }
  } else {
    for (var in variables) {
      tableName <- paste0("stem_and_leaf_", var)
      if (is.null(container[[tableName]])) {
        container[[tableName]] <- .descriptivesStemAndLeafCreateSingleTable(dataset[[var]], var, scale, width, atom)
        container[[tableName]]$dependOn(optionContainsValue = list(variables = var))
      }
    }
  }
}

.descriptivesStemAndLeafCreateSingleTable <- function(x, title, scale = 1, width = 80, atom = 1e-08) {
  tb <- createJaspTable(title = title)
  tb$addColumnInfo(name = "left", title = gettext("Stem"), type = "integer")
  tb$addColumnInfo(name = "sep", title = "", type = "separator")
  tb$addColumnInfo(name = "right", title = gettext("Leaf"), type = "string")

  if (length(na.omit(x)) < 2L) {
    tb$setError(gettext("A stem and leaf table could not be made because there are fewer than 2 non-missing observations"))
    return(tb)
  }

  # NOTE: graphics::stem is fast because it works in C, but it prints directly to the R output and returns NULL...
  # so we resort to capturing the string and manipulating it.
  # as.numeric ensures factors are handled correctly
  temp <- capture.output(graphics::stem(x, scale, width, atom))
  other <- temp[4:max(4, (length(temp) - 1L))]

  # parse the footnote so that we can translate it.
  # see https://github.com/wch/r-source/blob/c6710da21a5d869cd02889352483919ec3487000/src/library/graphics/src/stem.c#L105-L110 for details
  originalFootnote <- temp[2L]
  if (grepl("at the |", originalFootnote, fixed = TRUE)) {
    footnote <- gettext("The decimal point is at the |")
  } else {
    digits <- as.numeric(regmatches(originalFootnote, gregexpr("[[:digit:]]+", originalFootnote))[[1L]])

    footnote <- if (grepl("right", originalFootnote)) {
      sprintf(ngettext(
        digits,
        "The decimal point is %s digit to the right of the |",
        "The decimal point is %s digits to the right of the |"
      ), digits)
    } else {
      sprintf(ngettext(
        digits,
        "The decimal point is %s digit to the left of the |",
        "The decimal point is %s digits to the left of the |"
      ), digits)
    }
  }

  text <- strsplit(other, " | ", fixed = TRUE)
  left <- vapply(text, `[`, 1L, FUN.VALUE = character(1L))
  right <- vapply(text, function(x) if (length(x) == 1L) "" else x[2L], FUN.VALUE = character(1L))

  tb[["left"]] <- left
  tb[["sep"]] <- rep("|", length(left))
  tb[["right"]] <- right

  tb$addFootnote(footnote)

  return(tb)
}

.descriptivesHeatmaps <- function(container, dataset, variables, options) {
  axesNames <- c(options[["heatmapHorizontalAxis"]], options[["heatmapVerticalAxis"]])

  # we are not ready to plot
  if (length(variables) == 0 || isFALSE(options[["heatmapPlot"]]))
    return()

  axes <- dataset[axesNames]

  for (i in seq_along(variables)) {
    variableName <- variables[[i]]
    variableType <- options[["variables.types"]][[i]]
    variable <- dataset[[variableName]]

    if (options[["splitBy"]] == "") {
      .descriptivesCreateSingleHeatmap(container,
                                       axes,
                                       axesNames,
                                       variable,
                                       variableType,
                                       variableName,
                                       i, options)
    } else {
      container[[variableName]] <- createJaspContainer(variableName)
      splitBy <- dataset[, options[["splitBy"]]]
      groups <- levels(splitBy)

      for (g in seq_along(groups)) {
        activeCases <- groups[g] == splitBy
        .descriptivesCreateSingleHeatmap(container[[variableName]],
                                         axes[activeCases, ],
                                         axesNames,
                                         variable[activeCases],
                                         variableType,
                                         groups[g],
                                         g, options)
      }
    }
  }
}

.descriptivesCreateSingleHeatmap <- function(container, axes, axesNames, variable, variableType, plotName, position, options) {
  if (variableType == "scale") {
    data <- .descriptivesHeatmapAggregateData(variable, axes, options[["heatmapStatisticContinuous"]])
  } else {
    data <- .descriptivesHeatmapAggregateData(variable, axes, options[["heatmapStatisticDiscrete"]])
  }

  nLevels <- c(nlevels(data[["horizontal"]]), nlevels(data[["vertical"]]))
  plotSize <- c(200 + nLevels * 25) * c(options[["heatmapTileWidthHeightRatio"]], 1)
  if (any(plotSize > 700)) {
    plotSize <- (plotSize / max(plotSize)) * 700
  }
  if (options[["heatmapLegend"]]) {
    plotSize <- plotSize + c(100, 100)
  }

  if (any(table(data[["horizontal"]], data[["vertical"]]) > 1)) {
    jaspPlot <- createJaspPlot(title = plotName, error = gettext("There must be a unique value per combination of levels of horizontal and vertical axis!"), position = position)
  } else {
    jaspPlot <- createJaspPlot(title = plotName, width = plotSize[1], height = plotSize[2], position = position)

    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = horizontal, y = vertical, fill = value)) +
      ggplot2::geom_tile(color = "black", size = 1) +
      ggplot2::xlab(axesNames[1]) +
      ggplot2::ylab(axesNames[2]) +
      ggplot2::coord_fixed(ratio = 1 / options[["heatmapTileWidthHeightRatio"]])

    if (options[["heatmapDisplayValue"]]) {
      plot <- plot + ggplot2::geom_text(ggplot2::aes(x = horizontal, y = vertical, label = label),
        size = 8 * options[["heatmapDisplayValueRelativeTextSize"]]
      )
    }

    palette <- options[["colorPalette"]]
    # This factor check refers to the result of the summary statistic (mode/n/etc)
    plot <- plot + if (is.factor(data[["value"]])) {
      jaspGraphs::scale_JASPfill_discrete(palette)
    } else {
      jaspGraphs::scale_JASPfill_continuous(palette)
    }

    plot <- jaspGraphs::themeJasp(plot, legend.position = if (options[["heatmapLegend"]]) "right" else "none")

    jaspPlot$plotObject <- plot
  }

  container[[plotName]] <- jaspPlot
}

.descriptivesHeatmapAggregateData <- function(variable, groups, fun = c("identity", "meanArithmetic", "median", "mode", "length")) {
  fun <- match.arg(fun)
  mode <- function(x) {
    levels <- levels(x)
    out <- names(which.max(table(x)))
    factor(out, levels = levels)
  }
  data <- switch(fun,
    identity = cbind(groups, variable = variable),
    mean     = stats::aggregate(x = variable, by = groups, mean),
    median   = stats::aggregate(x = variable, by = groups, median),
    mode     = stats::aggregate(x = variable, by = groups, mode),
    length   = stats::aggregate(x = variable, by = groups, length)
  )

  colnames(data) <- c("horizontal", "vertical", "value")
  data$label <- if (is.numeric(data$value)) round(data$value, digits = 2) else data$value

  return(data)
}

.descriptivesLikertPlots <- function(dataset, name, options) {
  if(ncol(dataset) == 0) return()
  variables <- names(dataset)

  if (options[["likertPlotAssumeVariablesSameLevel"]]) {
    plotR <- .descriptivesLikertPlotsFill(dataset, variables, name, options)
    return(plotR)
  } else {
    plotResult <- createJaspContainer(title = name)
    for (var in variables) {
      data <- dataset[, names(dataset) %in% c(var), drop = FALSE]
      plotResult[[var]] <- .descriptivesLikertPlotsFill(data, var, var, options)
    }
    return(plotResult)
  }
}

.descriptivesLikertPlotsFill <- function(dataset, variables, name, options) {
  leng <- length(dataset)
  depends <- c("likertPlot", "splitBy", "variables", "likertPlotAssumeVariablesSameLevel", "likertPlotAdjustableFontSize")

  likPlot <- createJaspPlot(title = name, dependencies = depends, width = 1300, height = if (leng == 1) 250 else 200 * (leng * 0.8))
  errorMessage <- .descriptivesCheckPlotErrors(dataset, variables, obsAmount = "< 2")
  if (!is.null(errorMessage)) {
    likPlot$setError(gettextf("Plotting not possible: %s", errorMessage))
    return(likPlot)
  }

  # Likert Part: Preparing & summarize data in the likert format (% of levels per variable)
  nLevels <- nlevels(dataset[, 1])
  center <- (nLevels - 1) / 2 + 1
  lowRange <- 1:floor(center - 0.5)
  highRange <- ceiling(center + 0.5):nLevels

  if (!all(sapply(dataset, function(x) nlevels(x)) == nLevels)) {
    likPlot$setError(gettext("All categorical variables must have the same number of levels!"))
    return(likPlot)
  }
  if (center < 1.5) {
    likPlot$setError(gettext("Items must have 2 or more levels!"))
    return(likPlot)
  }

  # Summarizing item contribution and get clean data
  results <- data.frame()
  results <- data.frame(Response = seq_len(nLevels))
  for (i in seq_len(ncol(dataset))) {
    t <- table(dataset[, i])
    t <- (t / sum(t) * 100)
    results <- cbind(results, as.data.frame(t)[, 2])
    names(results)[ncol(results)] <- names(dataset)[i]
  }
  results <- as.data.frame(t(results))
  names(results) <- levels(dataset[, 1])
  results <- results[2:nrow(results), ]

  # Summarizing likert data: high, low, neutral %
  resultsTwo <- data.frame(
    Item = row.names(results),
    low = rep(NA, nrow(results)),
    neutral = rep(NA, nrow(results)),
    high = rep(NA, nrow(results))
  )
  resultsTwo$low <- if (length(lowRange) == 1) {
    results[, lowRange]
  } else {
    rowSums(results[, lowRange])
  }
  resultsTwo$high <- if (length(highRange) == 1) {
    results[, highRange]
  } else {
    rowSums(results[, highRange])
  }
  if (lowRange[length(lowRange)] + 1 != highRange[1]) {
    resultsTwo$neutral <- results[, (highRange[1] - 1)]
  }
  row.names(resultsTwo) <- seq_len(nrow(resultsTwo))

  results <- cbind(row.names(results), results)
  names(results)[1] <- "Item"
  row.names(results) <- seq_len(nrow(results))

  # Correcting for missing values in "results"
  for (i in 2:ncol(results)) {
    narows <- which(is.na(results[, i]))
    if (length(narows) > 0) {
      results[narows, i] <- 0
    }
  }
  # Correcting for missing values in "resultsTwo"
  narows <- which(is.na(resultsTwo$low))
  if (length(narows) > 0) {
    resultsTwo[narows, ]$low <- 0
  }
  narows <- which(is.na(resultsTwo$neutral))
  if (length(narows) > 0) {
    resultsTwo[narows, ]$neutral <- 0
  }
  narows <- which(is.na(resultsTwo$high))
  if (length(narows) > 0) {
    resultsTwo[narows, ]$high <- 0
  }

  lik <- list(results = results, levels = levels(dataset[, 1]), sum = resultsTwo)

  # Likert Plot Part:
  textSize <- 6
  textColor <- "black"
  yMin <- -100
  yMax <- 100
  yBuffer <- 5

  palette <- c(
    "#D8B365", "#E1C58B", "#EBD9B2", "#F5ECD8",
    "#D5ECEA", "#ACD9D5", "#83C6C0", "#5AB4AC"
  )
  cols <- scales::gradient_n_pal(palette, values = NULL)(seq(0, 1, length.out = nLevels))
  if (center %% 1 == 0) {
    cols[center] <- "grey90"
  }

  resultsLong <- stats::reshape(
    data = lik$results,
    idvar = "Item",
    v.name = c("value"),
    varying = c(names(lik$results[, 2:length(lik$results)])),
    times = c(names(lik$results[, 2:length(lik$results)])),
    timevar = "variable",
    new.row.names = seq_len(length(lik$results[2:length(lik$results)]) * length(lik$results$Item)),
    direction = "long"
  )

  resultsLong$Item <- factor(resultsLong$Item, levels = rev(lik$results$Item))
  orderLeg <- lik$levels # important for the correct legend sequence
  resultsLong$variable <- factor(resultsLong$variable, levels = orderLeg)

  # Make high, low, neutral values distinguishable
  rows <- which(resultsLong$variable %in% names(lik$results)[2:(length(lowRange) + 1)])
  resultsLong[rows, "value"] <- -1 * resultsLong[rows, "value"]
  if (center %% 1 == 0) {
    rowsMid <- which(resultsLong$variable %in% names(lik$results)[center + 1])
    tmp <- resultsLong[rowsMid, ]
    tmp$value <- tmp$value / 2 * -1
    resultsLong[rowsMid, "value"] <- resultsLong[rowsMid, "value"] / 2
    resultsLong <- rbind(resultsLong, tmp)
  }
  resultsLow <- resultsLong[resultsLong$value < 0, ]
  resultsHigh <- resultsLong[resultsLong$value > 0, ]

  p <- ggplot2::ggplot(resultsLong, ggplot2::aes(y = value, x = Item, group = Item)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey90") +
    ggplot2::geom_bar(data = resultsLow[nrow(resultsLow):1, ], ggplot2::aes(fill = variable), stat = "identity") +
    ggplot2::geom_bar(data = resultsHigh, ggplot2::aes(fill = variable), stat = "identity")

  names(cols) <- levels(resultsLong$variable)
  p <- p + ggplot2::scale_fill_manual("Response", breaks = names(cols), values = cols, drop = FALSE)

  p <- p + ggplot2::geom_text(
    data = lik$sum, # plot percent low
    y = yMin,
    ggplot2::aes(x = Item, label = paste0(round(low), "%")),
    size = textSize,
    hjust = 0.7,
    color = textColor
  ) +
    ggplot2::geom_text(
      data = lik$sum, # plot percent high
      y = 100,
      ggplot2::aes(x = Item, label = paste0(round(high), "%")),
      size = textSize,
      hjust = 0.3,
      color = textColor
    )

  if (nLevels %% 2 == 1) { # plot percent neutral
    p <- p + ggplot2::geom_text(
      data = lik$sum,
      y = 0,
      ggplot2::aes(x = Item, label = paste0(round(neutral), "%")),
      size = textSize,
      hjust = 0.5,
      color = textColor
    )
  }
  p <- p + ggplot2::coord_flip() +
    ggplot2::ylab("Percentage") +
    ggplot2::xlab("") +
    ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = function(x) {
      return(abs(x))
    }, limits = c(yMin - yBuffer, yMax + yBuffer)) +
    ggplot2::theme(legend.position = "bottom") + ggplot2::guides(fill = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::theme(panel.background = ggplot2::element_rect(size = 1, color = "grey90", fill = NA)) +
    ggplot2::theme(text = ggplot2::element_text(size = 22.5), axis.title.x = ggplot2::element_text(size = 18))

  p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(
    size = switch(options[["likertPlotAdjustableFontSize"]],
      "small"  = 20,
      "medium" = 22.5,
      "large"  = 25
    )
  ))

  likPlot$plotObject <- p

  return(likPlot)
}

.descriptivesParetoPlots <- function(dataset, variable, options) {
  if (options[["splitBy"]] != "") {
    split <- names(dataset)

    plotResult <- createJaspContainer(title = variable)
    plotResult$dependOn(options = "splitBy", optionContainsValue = list(variables = variable))

    for (i in split) {
      plotResult[[i]] <- .descriptivesParetoPlots_SubFunc(dataset[[i]], variable, i, options)
      plotResult[[i]]$dependOn(optionsFromObject = plotResult)
    }
    return(plotResult)

  } else {
    pPlot <- .descriptivesParetoPlots_SubFunc(dataset, variable, variable, options)
    pPlot$dependOn(options = "splitBy", optionContainsValue = list(variables = variable))

    return(pPlot)
  }
}

.descriptivesParetoPlots_SubFunc <- function(dataset, variable, title, options) {
  parPlot <- createJaspPlot(title = title, width = options[["plotWidth"]], height = options[["plotHeight"]])

  countVar <- if (options[["paretoAddCountVariable"]] != "") options[["paretoAddCountVariable"]] else NULL

  errorMessage <- .descriptivesCheckPlotErrors(dataset, c(variable, countVar),
                                                    obsAmount = "< 2")

  if (options[["paretoAddCountVariable"]] == "") {
    # single column as vector
    column <- dataset[[variable]]
    column <- column[!is.na(column)]
  } else { # there is a count variable
    # check errors for the count variable
    if (variable == options[["paretoAddCountVariable"]]) {
      parPlot$setError(gettext("The variable used for counts must be different from the variable used for categories."))
      return(parPlot)
    }
    column <- dataset[c(variable, options[["paretoAddCountVariable"]])]
    column <- column[stats::complete.cases(column), , drop = FALSE]
  }
  if (!is.null(errorMessage)) {
    parPlot$setError(gettextf("Plotting not possible: %s", errorMessage))
  } else if (length(column) > 0) {
    parPlot$plotObject <- .descriptivesParetoPlots_Fill(column, variable, options)
  }

  return(parPlot)
}

.descriptivesParetoPlots_Fill <- function(column, variable, options) {

  # tableNote <- ""
  # Normalize input: either a vector OR a 2-col data.frame (category, COUNT)
  if (is.data.frame(column) && ncol(column) == 2) {
    # does the x-variable have duplicated categories, collapse them
    if (length(unique(column[[1]])) != nrow(column)) {
      column <- aggregate(as.numeric(as.character(column[[2]])), by = list(column[[1]]), FUN = sum)
      # tableNote <- gettext("Note: Categories with multiple counts have been collapsed.")
    }
    tb <- data.frame(level = column[[1]], Freq = as.numeric(as.character(column[[2]])))
  } else {
    tab <- table(column, useNA = "no")
    tb  <- data.frame(level = names(tab), Freq = as.numeric(tab))
  }

  # Order & cumulative %
  tb <- tb[order(tb$Freq, decreasing = TRUE), , drop = FALSE]
  rownames(tb) <- NULL
  tb$cums <- cumsum(tb$Freq)
  tb$cums <- 100 * tb$cums / tail(tb$cums, 1)

  # Shared quantities
  yBreaks    <- jaspGraphs::getPrettyAxisBreaks(c(0, tb$Freq))
  scaleRight <- max(tb$cums) / max(yBreaks)

  # --- Unified Pareto line + guides (left axis = % when shifted) ---
  shift    <- options[["paretoShiftAccumulationLine"]]
  tot      <- sum(tb$Freq)
  cumC     <- cumsum(tb$Freq)
  cumP     <- cumC / tot * 100
  xFac     <- factor(tb$level, levels = tb$level)
  barWidth <- 0.9
  edgeShift <- barWidth / 2

  # After you build tb (already sorted by Freq desc):
  tb$level <- factor(tb$level, levels = tb$level)  # lock global order

  yLabel <- gettext("Counts")
  yBar <- tb$Freq
  p <- ggplot2::ggplot(tb, ggplot2::aes(x = level, y = yBar)) +
    ggplot2::geom_col(fill = "grey", col = "black", width = 0.9) +
    ggplot2::xlab(variable) +
    ggplot2::ylab(yLabel) +
    ggplot2::scale_x_discrete(limits = levels(tb$level), drop = FALSE)

  # Line: cumulative counts when shifted, otherwise the old scaled % line
  yLine <- if (shift) cumC else tb$cums / scaleRight
  dfL   <- data.frame(level = tb$level, y = yLine)

  # left-axis breaks
  yBreaksLeft <- if (shift) {
    br <- jaspGraphs::getPrettyAxisBreaks(c(0, tot))
    # drop any break too close to tot to avoid overlapping tick labels
    step <- if (length(br) > 1) min(diff(br)) else tot
    br <- br[abs(br - tot) > step * 0.5 | br == 0]
    sort(unique(c(br, tot)))
  } else {
    yBreaks
  }

  p <- p +
    ggplot2::geom_line(
      data = dfL,
      ggplot2::aes(x = level, y = y, group = 1),
      colour = "black", linewidth = 1, inherit.aes = FALSE,
      position = if (shift) ggplot2::position_nudge(x = edgeShift) else "identity"
    ) +
    ggplot2::geom_point(
      data = dfL,
      ggplot2::aes(x = level, y = y),
      colour = "steelblue", size = 3, inherit.aes = FALSE,
      position = if (shift) ggplot2::position_nudge(x = edgeShift) else "identity"
    ) +
    ggplot2::scale_y_continuous(
      name   = yLabel,
      breaks = yBreaksLeft,
      limits = if (shift) c(0, tot) else range(yBreaksLeft),
      sec.axis = ggplot2::sec_axis(
        transform = if (shift) ~ . * 100 / tot else ~ . * scaleRight,
        name      = gettextf("Percentage (%%)"),
        breaks    = seq(0, 100, 20)
      )
    ) +
    jaspGraphs::geom_rangeframe(sides = "rbl") +
    jaspGraphs::themeJaspRaw() +
    ggplot2::theme(plot.margin = ggplot2::margin(5))

  # Pareto-rule guides (e.g., 95%) – vertical up from x-axis, then toward the % axis
  if (options[["paretoPlotRule"]]) {
    perc   <- options[["paretoPlotRuleCi"]]   # e.g., 0.95
    target <- perc * 100                      # % value

    # find crossing on the cumulative % curve
    k <- which(cumP >= target)[1]
    xPos <- if (is.na(k) || k <= 1) {
      if (shift) edgeShift else 1
    } else {
      p0 <- cumP[k - 1]; p1 <- cumP[k]
      base <- (k - 1) + (target - p0) / (p1 - p0)
      if (shift) base + edgeShift else base
    }

    # y height in *left-axis* units (counts)
    # Map percentage to count scale
    yGuideLeft <- if (shift) perc * tot else perc * max(yBreaks)

    # horizontal endpoint: always to the right where percentage axis lives
    xEnd <- length(levels(tb$level)) + 0.5

    p <- p +
      # vertical from x-axis up to the % height
      ggplot2::geom_segment(
        ggplot2::aes(x = xPos, xend = xPos, y = 0, yend = yGuideLeft),
        linetype = "dashed", color = "orange", linewidth = 1.3, inherit.aes = FALSE
      ) +
      # horizontal from the crossing toward the % axis side
      ggplot2::geom_segment(
        ggplot2::aes(x = xPos, xend = xEnd, y = yGuideLeft, yend = yGuideLeft),
        linetype = "dashed", color = "orange", linewidth = 1.3, inherit.aes = FALSE
      )
  }

  # possibly tilt the x-axis labels
  if (options[["paretoPlotTiltXAxisLabels"]]) {
    levs <- decodeColNames(levels(tb$level))
    maxChars <- max(nchar(levs))
    titlePad <- max(10, min(40, ceiling(maxChars * 0.6)))  # space to the title
    tickPad  <- 8  # space from ticks to labels (pts) — bump to taste

    p <- p +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(
          angle  = 45, hjust = 1, vjust = 1,
          margin = ggplot2::margin(t = tickPad)   # <-- key line
        ),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = titlePad)),
        plot.margin  = ggplot2::margin(t = 5, r = 5, b = titlePad + tickPad + 6, l = 5)
      ) +
      ggplot2::coord_cartesian(clip = "off")
  }

  return(p)


}
.descriptivesDensityPlots <- function(container, dataset, variables, options) {
  # we are not ready to plot
  if (length(variables) == 0)
    return()

  if (options[["densityPlotSeparate"]] != "")
      separator <- as.factor(dataset[[options[["densityPlotSeparate"]]]])

  for (i in seq_along(variables)) {

    variableName <- variables[[i]]
    variable <- dataset[[variableName]]
    variableType <- options[["variables.types"]][[i]]
    data <- if (options[["densityPlotSeparate"]] != "") {
      data.frame(variable, separator)
    } else {
      data.frame(variable)
    }

    if (options[["splitBy"]] != "") {
      container[[variableName]] <- createJaspContainer(variableName)
      data$split <-  dataset[[options[["splitBy"]]]]
      data <- na.omit(data)
      groups <- levels(data$split)

      for (g in seq_along(groups)) {
        active <- groups[g] == data$split
        .descriptivesDensityPlotsFill(container[[variableName]], data[active, ], groups[g], variableName, g, options, variableType)
      }
    } else {
      data <- na.omit(data)
      .descriptivesDensityPlotsFill(container, data, variableName, variableName, i, options, variableType)
    }
  }
}

.descriptivesDensityPlotsFill <- function(container, data, plotName, axeName, position, options, variableType) {
  jaspGraphs::graphOptions(palette = options[["colorPalette"]])
  densPlot <- createJaspPlot(title = plotName, width = 480, height = 320, position = position)

  errorMessage <- .descriptivesCheckPlotErrors(data, colnames(data), obsAmount = "< 2")
  if (!is.null(errorMessage)) {

      densPlot$setError(errorMessage)

  } else if (options[["densityPlotSeparate"]] != "" && any(table(data$separator) == 1)) {

    densPlot$setError(gettext("Levels within variable require at least two or more data points!"))

  } else if (variableType == "nominal" || variableType == "ordinal")  {

    if (options[["densityPlotCategoricalType"]] == "condProp" && options[["densityPlotSeparate"]] != "") {
      tb <- as.data.frame(table(data) / rowSums(table(data), na.rm = TRUE))
      yAxeName <- "Conditional proportion"
    } else if (options[["densityPlotCategoricalType"]] == "count") {
      tb <- as.data.frame(table(data))
      yAxeName <- "Counts"
    } else { # if no separator but request cond prop, then give proportions
      tb <- as.data.frame(table(data) / sum(table(data), na.rm = TRUE))
      yAxeName <- "Proportion"
    }

    p <- if (options[["densityPlotSeparate"]] != "") {
      ggplot2::ggplot(data = tb, ggplot2::aes(x = variable, y = Freq, fill = separator))
    } else {
      ggplot2::ggplot(data = tb, ggplot2::aes(x = variable, y = Freq))
    }

    densPlot$plotObject <- p +
      ggplot2::geom_bar(stat = "identity", col = "black", width = 0.7, size = .3, position = ggplot2::position_dodge(width = 0.7)) +
      ggplot2::xlab(axeName) +
      ggplot2::scale_fill_manual(values = jaspGraphs::JASPcolors(options[["colorPalette"]]), name = options[["densityPlotSeparate"]]) +
      ggplot2::ylab(yAxeName) +
      jaspGraphs::geom_rangeframe() +
      jaspGraphs::themeJaspRaw() +
      ggplot2::theme(plot.margin =  ggplot2::margin(5), legend.position = "right")

  } else if (options[["densityPlotType"]] == "density") {

    densPlot$plotObject <- jaspGraphs::jaspHistogram(x = data[["variable"]],
                                                     histogram = FALSE,
                                                     density = TRUE,
                                                     densityColor = TRUE,
                                                     densityShade = TRUE,
                                                     densityShadeAlpha = 1 - (options[["densityPlotTransparency"]]/100),
                                                     xName = axeName,
                                                     groupingVariableName = options[["densityPlotSeparate"]],
                                                     groupingVariable = data[["separator"]],
                                                     histogramPosition = options[["customHistogramPosition"]])


  } else if (options[["densityPlotType"]] == "histogram") {

    densPlot$plotObject <- jaspGraphs::jaspHistogram(x = data[["variable"]],
                                                     xName = axeName,
                                                     groupingVariableName = options[["densityPlotSeparate"]],
                                                     groupingVariable = data[["separator"]],
                                                     binWidthType = "sturges",
                                                     histogramPosition = options[["customHistogramPosition"]])
    if (options[["customHistogramPosition"]] == "identity")
      densPlot$plotObject$layers[[1]]$aes_params$alpha <- 1 - (options[["densityPlotTransparency"]]/100)
  }


  container[[plotName]] <- densPlot
}

## Variable summary plots: nominal -------------------------------------------------------------------------------------

.descriptivesSummaryNominal <- function(container, dataset, options) {
  variables <- unlist(options[["summaryPlotNominalVariablesSel"]])

  # not ready to plot
  if (length(variables) == 0)
    return()

  for (i in seq_along(variables)) {
    variableName <- variables[[i]]

    # keep the NAs, they are needed for the Missing count
    variable <- dataset[[variableName]]
    freqTable <- .descriptivesSummaryNominalFreqTable(variable)

    if (options[["summaryNominalType"]] == "barPlot"){
      # many levels: the least frequent ones are merged into "Other" (the panel still describes all of them)
      plotTable <- .descriptivesSummaryNominalMergeLevels(freqTable, .descriptivesSummaryNominalMaxLevels - 1)
      height    <- if (.descriptivesSummaryNominalBarSlanted(plotTable, options)) 400 else 320
      plot <- createJaspPlot(title = gettextf("Summary for %s", variableName), width = 720, height = height, position = i)
      plot$dependOn(optionContainsValue = list(summaryPlotNominalVariablesSel = variableName))
      container[[variableName]] <- plot

      errorMessage <- .descriptivesCheckPlotErrors(dataset, variableName, obsAmount = "< 1")

      if (!is.null(errorMessage)){
        plot$setError(errorMessage)
        next
      }
      main <- .descriptivesSummaryNominalBarPlot(plotTable, options)
    }
      
    
    if (options[["summaryNominalType"]] == "stackedPlot") {
      plotTable <- .descriptivesSummaryNominalStackedMergeLevels(freqTable, options)
      # taller when some labels have to go above the bar
      needsRoom <- !all(.descriptivesSummaryNominalStackedLabels(plotTable, options)$inside)
      plot <- createJaspPlot(title = variableName, width = 720, height = if (needsRoom) 400 else 320, position = i)
      plot$dependOn(optionContainsValue = list(summaryPlotNominalVariablesSel = variableName))
      container[[variableName]] <- plot

      errorMessage <- .descriptivesCheckPlotErrors(dataset, variableName, obsAmount = "< 1")

      if (!is.null(errorMessage)){
        plot$setError(errorMessage)
        next
      }
      main <- .descriptivesSummaryNominalStackedPlot(plotTable, options)
    }

    nMerged <- nrow(freqTable) - nrow(plotTable) + 1
    if (nMerged > 1)
      main <- main +
        ggplot2::labs(caption = gettextf("Other: the %d levels with the lowest frequencies.", nMerged)) +
        ggplot2::theme(plot.caption = ggplot2::element_text(size = 12, hjust = 0))

    if (options[["summaryNominalDisp"]]) {
      panel <- .descriptivesSummaryNominalPanel(variable, freqTable)
      # margins between plot and panel, so it is clear which tick labels belong to the plot
      main  <- patchwork::wrap_plots(main + ggplot2::theme(plot.margin = ggplot2::margin(r = 12)),
                                     panel + ggplot2::theme(plot.margin = ggplot2::margin(l = 12)),
                                     widths = c(3, 1))
    }
    plot$plotObject <- main

  }
}


.descriptivesSummaryNominalFreqTable <- function(variable) {
  freqNominal <- table(variable)
  freqNominal <- sort(freqNominal, decreasing = TRUE)
  freqNominal <- freqNominal[as.vector(freqNominal) != 0]

  naCount <- sum(is.na(variable))
  if (naCount > 0)
    freqNominal[["Missing"]] <- naCount
  
  propNominal <- as.vector(freqNominal) / sum(as.vector(freqNominal))
  freqNominalFrame <- data.frame(level = names(freqNominal), 
                                 count = as.vector(unname(freqNominal)),
                                 prop = as.vector(propNominal))

  factor(freqNominalFrame$level, levels = freqNominalFrame$level)
  
  return(freqNominalFrame)
}

# Maximum number of bars, "Other" included, Missing not
.descriptivesSummaryNominalMaxLevels <- 8

# Keep the `keep` most frequent levels and sum the rest into "Other" (only if that merges at least two levels)
.descriptivesSummaryNominalMergeLevels <- function(freqTable, keep) {
  isMissing <- freqTable$level == "Missing"
  real      <- freqTable[!isMissing, ]
  if (nrow(real) <= keep + 1)
    return(freqTable)

  rest  <- real[-seq_len(keep), ]
  other <- data.frame(level = "Other", count = sum(rest$count), prop = sum(rest$prop))

  return(rbind(real[seq_len(keep), ], other, freqTable[isMissing, ]))
}

# Stacked plot has no legend: merge levels until every remaining segment can be labelled
.descriptivesSummaryNominalStackedMergeLevels <- function(freqTable, options) {
  barPx <- .descriptivesSummaryNominalStackedBarPx(options)
  for (keep in (.descriptivesSummaryNominalMaxLevels - 1):1) {
    merged  <- .descriptivesSummaryNominalMergeLevels(freqTable, keep)
    labels  <- .descriptivesSummaryNominalStackedLabels(merged, options)
    outside <- .descriptivesSummaryNominalOutsideLabels(labels[!labels$inside, ], barPx)
    if (nrow(outside) == sum(!labels$inside))
      break
  }

  return(merged)
}

# Bar plot with the count above each bar
.descriptivesSummaryNominalBarPlot <- function(freqTable, options) {
  freqTable$level <- factor(freqTable$level, levels = freqTable$level)
  p <- ggplot2::ggplot(freqTable, ggplot2::aes(x = level, y = count, fill = level)) +
       ggplot2::geom_bar(stat = "identity", col = "black", width = 0.7, linewidth = .3) +
       ggplot2::geom_label(ggplot2::aes(label = count), vjust = -0.3, size = .descriptivesSummaryNominalTextSize,
                           fill = "white", colour = "black", linewidth = .3, show.legend = FALSE)
  
  values <- .descriptivesSummaryNominalColors(freqTable, options)

  p <- p + 
       ggplot2::scale_fill_manual(values = values) + 
       ggplot2::theme(legend.position = "none") +
       jaspGraphs::geom_rangeframe() + 
       jaspGraphs::themeJaspRaw() +
       ggplot2::ylab(gettext("Frequency")) +
       ggplot2::xlab(gettext("Level"))
  
  if (.descriptivesSummaryNominalBarSlanted(freqTable, options))
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))

  # headroom for the count label above the tallest bar
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(freqTable$count) * 1.35))
  p <- p + 
       ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks))

  return(p)
}

# Slant the level names when the longest does not fit under its bar (~11 px per character)
.descriptivesSummaryNominalBarSlanted <- function(freqTable, options) {
  slotPx <- ((if (options[["summaryNominalDisp"]]) 720 * 3 / 4 else 720) - 100) / nrow(freqTable)
  return(max(nchar(as.character(freqTable$level))) * 11 > slotPx)
}

# Stacked plot: one horizontal 100% bar, each segment labelled with level, count and percentage
.descriptivesSummaryNominalStackedPlot <- function(freqTable, options) {
  freqTable$level <- factor(freqTable$level, levels = freqTable$level)
  values <- .descriptivesSummaryNominalColors(freqTable, options)

  textSize  <- .descriptivesSummaryNominalTextSize
  barPx     <- .descriptivesSummaryNominalStackedBarPx(options)
  freqTable <- .descriptivesSummaryNominalStackedLabels(freqTable, options)

  inside  <- freqTable[freqTable$inside, ]
  outside <- .descriptivesSummaryNominalOutsideLabels(freqTable[!freqTable$inside, ], barPx)

  p <- ggplot2::ggplot(freqTable, ggplot2::aes(x = prop, y = 1, fill = level)) +
       # orientation = "y" is needed because y is numeric
       ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE), orientation = "y",
                    col = "black", width = 0.5, linewidth = .3) +
       .descriptivesSummaryNominalStackedBoxes(inside, ggplot2::aes(x = centre, y = 1), textSize)

  # labels that do not fit inside their segment go above the bar, with a line to the segment
  if (nrow(outside) > 0)
    p <- p +
      ggplot2::geom_segment(data = outside, ggplot2::aes(x = centre, xend = labelX, y = 1.26, yend = 1.4),
                            colour = "grey50", linewidth = .3, inherit.aes = FALSE) +
      .descriptivesSummaryNominalStackedBoxes(outside, ggplot2::aes(x = labelX, y = 1.4), textSize, vjust = 0)

  p <- p +
       ggplot2::scale_fill_manual(values = values, name = NULL) +
       # limits in coord_cartesian(): the proportions can sum to slightly over 1, which would drop the last segment
       ggplot2::scale_x_continuous(breaks = seq(0, 1, .25), labels = function(x) paste0(x * 100, "%")) +
       ggplot2::xlab(gettext("Percentage")) + ggplot2::ylab(NULL) +
       jaspGraphs::geom_rangeframe(sides = "b") +
       ggplot2::scale_y_continuous(limits = c(0.7, if (nrow(outside) > 0) 2 else 1.3), expand = c(0, 0)) +
       ggplot2::coord_cartesian(xlim = c(0, 1), clip = "off") +
       jaspGraphs::themeJaspRaw(legend.position = "none") +
       ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())
  
  return(p)
}

# Boxed labels with a bold level name: one label cannot mix font faces, so the box and the two parts are layered
.descriptivesSummaryNominalStackedBoxes <- function(data, mapping, textSize, vjust = 0.5) {
  list(
    ggplot2::geom_label(data = data, mapping = utils::modifyList(mapping, ggplot2::aes(label = label)),
                        size = textSize, fontface = "bold", fill = "white", text.colour = NA,
                        border.colour = "black", linewidth = .3, vjust = vjust, lineheight = .9, inherit.aes = FALSE),
    ggplot2::geom_label(data = data, mapping = utils::modifyList(mapping, ggplot2::aes(label = labelName)),
                        size = textSize, fontface = "bold", fill = NA, colour = "black", border.colour = NA,
                        vjust = vjust, lineheight = .9, inherit.aes = FALSE),
    ggplot2::geom_label(data = data, mapping = utils::modifyList(mapping, ggplot2::aes(label = labelNumbers)),
                        size = textSize, fill = NA, colour = "black", border.colour = NA,
                        vjust = vjust, lineheight = .9, inherit.aes = FALSE)
  )
}

# Text size of the boxed numbers, and the stacked bar's width in px (measured)
.descriptivesSummaryNominalTextSize <- 6.5
.descriptivesSummaryNominalStackedBarPx <- function(options) {
  if (options[["summaryNominalDisp"]]) 720 * 3 / 4 - 85 else 720 - 60
}

# Labels of the stacked plot, each segment's centre (0-1) and whether its label fits inside it
.descriptivesSummaryNominalStackedLabels <- function(freqTable, options) {
  textSize <- .descriptivesSummaryNominalTextSize
  freqTable$label        <- sprintf("%s\n%d\n%.1f%%", freqTable$level, freqTable$count, freqTable$prop * 100)
  freqTable$labelName    <- sprintf("%s\n \n ", freqTable$level)
  freqTable$labelNumbers <- sprintf(" \n%d\n%.1f%%", freqTable$count, freqTable$prop * 100)
  freqTable$widthPx <- pmax(nchar(as.character(freqTable$level)) * 1.1, nchar(freqTable$count),
                            nchar(sprintf("%.1f%%", freqTable$prop * 100))) *
    textSize * 0.6 * 96 / 25.4 + 16
  freqTable$centre  <- cumsum(freqTable$prop) - freqTable$prop / 2
  freqTable$inside  <- freqTable$prop * .descriptivesSummaryNominalStackedBarPx(options) >= freqTable$widthPx

  return(freqTable)
}

# Positions (labelX) of the labels above the bar: as close to their segment as possible without overlapping
# (isotonic regression = least total displacement). The smallest are dropped if they cannot all fit.
.descriptivesSummaryNominalOutsideLabels <- function(labels, barPx) {
  if (nrow(labels) == 0)
    return(cbind(labels, labelX = numeric(0)))

  keep <- order(labels$prop, decreasing = TRUE)
  keep <- keep[cumsum(labels$widthPx[keep]) <= barPx]
  labels <- labels[sort(keep), ]
  if (nrow(labels) == 0)
    return(cbind(labels, labelX = numeric(0)))

  w      <- labels$widthPx / barPx
  c      <- c(0, cumsum((w[-1] + w[-length(w)]) / 2))
  z      <- if (nrow(labels) == 1) labels$centre else stats::isoreg(seq_along(c), labels$centre - c)$yf
  # keep the row inside the bar at both ends
  z      <- pmin(pmax(z, w[1] / 2), max(w[1] / 2, 1 - c[length(c)] - w[length(w)] / 2))
  labels$labelX <- z + c

  return(labels)
}



.descriptivesSummaryNominalColors <- function(freqTable, options){
  realLevels <- setdiff(levels(freqTable$level), c("Other", "Missing"))
  palette <- jaspGraphs::JASPcolors(options[["colorPalette"]])
  # continuous palettes (viridis, inferno) are sampled evenly, discrete ones are used as they are
  cols <- if (length(palette) > 20) {
    jaspGraphs::JASPcolors(options[["colorPalette"]], asFunction = TRUE)(length(realLevels))
  } else {
    rep_len(palette, length(realLevels))
  }
  values <- stats::setNames(cols, realLevels)

  # light grey for Other, darker grey for Missing
  if ("Other" %in% levels(freqTable$level))
    values <- c(values, Other = "grey90")
  if ("Missing" %in% levels(freqTable$level))
    values <- c(values, Missing = "grey")
  
  return(values)
}

.descriptivesSummaryNominalPanel <- function(variable, freqTable) {
  nTotal <- length(variable)
  nMissing <- sum(is.na(variable))
  nValid <- nTotal - nMissing
  validTab <- freqTable[freqTable$level != "Missing",]

  rows <- data.frame(
    group = c(1, 1, 1, 2, 3, 3, 3),
    label = c(gettext("Total N"), gettext("Valid N"), gettext("Missing N"),
              gettext("Valid %"),
              gettext("Levels"), gettext("Mode"), gettext("Mode count")),
    value = c(nTotal, nValid, nMissing,
              sprintf("%.1f%%", nValid / nTotal * 100),
              nrow(validTab), as.character(validTab$level[1]), validTab$count[1]),
    stringsAsFactors = FALSE
  )

  rows$y <- -(seq_len(nrow(rows)) + (rows$group - 1) * 0.6)

  pPanel <- ggplot2::ggplot() +
    ggplot2::geom_text(data = rows, ggplot2::aes(x = 0, y = y, label = label), hjust = 0, size = 3.5, fontface = "bold") +
    ggplot2::geom_text(data = rows, ggplot2::aes(x = 1, y = y, label = value), hjust = 1, size = 3.5)

  seps <- tapply(rows$y, rows$group, min)[-3] - 0.6

  pPanel <- pPanel + 
            ggplot2::geom_hline(yintercept = seps, linewidth = .3) +
            ggplot2::scale_x_continuous(limits = c(-0.1, 1.1)) +
            ggplot2::scale_y_continuous(limits = c(min(rows$y) - 0.7, max(rows$y) + 0.7)) +
            ggplot2::theme_void() +
            ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .3))
  
  return(pPanel)
}

## Variable summary plots: ordinal -------------------------------------------------------------------------------------

.descriptivesSummaryOrdinal <- function(container, dataset, options) {
  variables <- unlist(options[["summaryPlotOrdinalVariablesSel"]])

  # not ready to plot
  if (length(variables) == 0)
    return()


  for (i in seq_along(variables)) {
    variableName <- variables[[i]]


    variable <- dataset[[variableName]]
    freqTable <- .descriptivesSummaryOrdinalFreqTable(variable)

    plot <- createJaspPlot(title = gettextf("Summary for %s", variableName), width = 720, height = 320, position = i)
    plot$dependOn(optionContainsValue = list(summaryPlotOrdinalVariablesSel = variableName))
    container[[variableName]] <- plot

    errorMessage <- .descriptivesCheckPlotErrors(dataset, variableName, obsAmount = "< 1")

    if (!is.null(errorMessage)){
      plot$setError(errorMessage)
      next
    }
    main <- .descriptivesSummaryOrdinalBarPlot(freqTable, options)
      

    if (options[["summaryOrdinalDisp"]]) {
      panel <- .descriptivesSummaryOrdinalPanel(variable, freqTable, options)
      main  <- patchwork::wrap_plots(main + ggplot2::theme(plot.margin = ggplot2::margin(r = 12)),
                                     panel + ggplot2::theme(plot.margin = ggplot2::margin(l = 12)),
                                     widths = c(3, 1))
    }
    plot$plotObject <- main

  }
}


.descriptivesSummaryOrdinalFreqTable <- function(variable) {
  freqOrdinal <- table(variable)

  naCount <- sum(is.na(variable))
  if (naCount > 0)
    freqOrdinal[["Missing"]] <- naCount
  
  propOrdinal <- as.vector(freqOrdinal) / sum(as.vector(freqOrdinal))
  freqOrdinalFrame <- data.frame(level = names(freqOrdinal), 
                                 count = as.vector(unname(freqOrdinal)),
                                 prop = as.vector(propOrdinal))

  factor(freqOrdinalFrame$level, levels = freqOrdinalFrame$level)
  
  return(freqOrdinalFrame)
}

# Bar plot with the count above each bar
.descriptivesSummaryOrdinalBarPlot <- function(freqTable, options) {
  freqTable$level <- factor(freqTable$level, levels = freqTable$level)
  p <- ggplot2::ggplot(freqTable, ggplot2::aes(x = level, y = count, fill = level)) +
       ggplot2::geom_bar(stat = "identity", col = "black", width = 0.7, linewidth = .3) +
       ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.5, size = 5.5)

  values <- .descriptivesSummaryOrdinalColors(freqTable, options)

  p <- p + 
       ggplot2::scale_fill_manual(values = values) + 
       ggplot2::theme(legend.position = "none") +
       jaspGraphs::geom_rangeframe() + 
       jaspGraphs::themeJaspRaw() +
       ggplot2::ylab(gettext("Frequency")) +
       ggplot2::xlab(gettext("Level"))
  
  # headroom for the count label above the tallest bar
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(freqTable$count) * 1.15))
  p <- p + 
       ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks))

  return(p)
}



.descriptivesSummaryOrdinalColors <- function(freqTable, options) {
  levelNames <- as.character(freqTable$level)
  realLevels <- setdiff(levelNames, "Missing")
  nLevels    <- length(realLevels)

  base   <- jaspGraphs::JASPcolors(options[["colorPalette"]])[1]
  cols   <- grDevices::colorRampPalette(c("white", base))(nLevels + 2)[-(1:2)]
  values <- stats::setNames(cols, realLevels)

  if ("Missing" %in% levelNames)
    values <- c(values, Missing = "grey")

  return(values)
}

.descriptivesSummaryOrdinalPanel <- function(variable, freqTable, options) {
  nTotal <- length(variable)
  nMissing <- sum(is.na(variable))
  nValid <- nTotal - nMissing
  validTab <- freqTable[freqTable$level != "Missing", ]

  modeIndex <- which.max(validTab$count)

  qType       <- .ordinalQuantileType(as.integer(options[["quantilesType"]]))
  ordQuantile <- function(p) as.character(quantile(na.omit(variable), p, names = FALSE, type = qType))

  rows <- data.frame(
    group = c(1, 1, 1, 2, 3, 3, 3, 4, 4, 4),
    label = c(gettext("Total N"), gettext("Valid N"), gettext("Missing N"),
              gettext("Valid %"),
              gettext("Levels"), gettext("Mode"), gettext("Mode count"),
              gettext("Median"), gettext("25% quantile"), gettext("75% quantile")),

    value = c(nTotal, nValid, nMissing,
              sprintf("%.1f%%", nValid / nTotal * 100),
              nrow(validTab), as.character(validTab$level[modeIndex]), validTab$count[modeIndex],
              ordQuantile(0.5), ordQuantile(0.25), ordQuantile(0.75)),
    stringsAsFactors = FALSE
  )

  rows$y <- -(seq_len(nrow(rows)) + (rows$group - 1) * 0.6)

  pPanel <- ggplot2::ggplot() +
    ggplot2::geom_text(data = rows, ggplot2::aes(x = 0, y = y, label = label), hjust = 0, size = 3.5, fontface = "bold") +
    ggplot2::geom_text(data = rows, ggplot2::aes(x = 1, y = y, label = value), hjust = 1, size = 3.5)

  seps <- tapply(rows$y, rows$group, min)
  seps <- seps[-length(seps)] - 0.6

  pPanel <- pPanel +
            ggplot2::geom_hline(yintercept = seps, linewidth = .3) +
            ggplot2::scale_x_continuous(limits = c(-0.1, 1.1)) +
            ggplot2::scale_y_continuous(limits = c(min(rows$y) - 0.7, max(rows$y) + 0.7)) +
            ggplot2::theme_void() +
            ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .3))

  return(pPanel)
}


## Variable summary plots: scale ----------------------------------------------------------------------------------------
# Main plot (raincloud or histogram) plus optional panels: time series, Q-Q plot, numbers

.descriptivesSummaryScale <- function(container, dataset, options) {
  variables <- unlist(options[["summaryPlotScaleVariablesSel"]])

  # not ready to plot
  if (length(variables) == 0)
    return()

  # the order is the reading order of the layout
  showPanels <- c(timeSeries = options[["summaryScaleTimeSeries"]],
                  qq         = options[["summaryScaleQQ"]],
                  numbers    = options[["summaryScaleNumbers"]])
  size     <- .descriptivesSummaryScaleSize(showPanels)
  mainWidth <- .descriptivesSummaryScaleMainWidth(size[["width"]], names(showPanels)[showPanels])

  for (i in seq_along(variables)) {
    variableName <- variables[[i]]

    # keep the NAs, they are needed for the Missing count
    variable <- dataset[[variableName]]

    plot <- createJaspPlot(title = gettextf("Summary for %s", variableName),
                           width = size[["width"]], height = size[["height"]], position = i)
    plot$dependOn(optionContainsValue = list(summaryPlotScaleVariablesSel = variableName))
    container[[variableName]] <- plot

    # the density and the Q-Q plot need at least 2 observations
    errorMessage <- .descriptivesCheckPlotErrors(dataset, variableName, obsAmount = "< 2")

    if (!is.null(errorMessage)) {
      plot$setError(errorMessage)
      next
    }

    main <- if (options[["summaryScaleMainType"]] == "rain") {
      .descriptivesSummaryScaleRaincloud(variable, options, mainWidth)
    } else {
      .descriptivesSummaryScaleHistogram(variable, options, mainWidth)
    }

    panels <- list()
    if (showPanels[["timeSeries"]])
      panels[["timeSeries"]] <- .descriptivesSummaryScaleTimeSeries(variable, options, variableName,
                                                                     rowNumbers = as.integer(row.names(dataset)))
    if (showPanels[["qq"]])         panels[["qq"]]         <- .descriptivesSummaryScaleQQ(variable, options)
    if (showPanels[["numbers"]])    panels[["numbers"]]    <- .descriptivesSummaryScalePanel(variable, options)

    mainTitle <- if (options[["summaryScaleMainType"]] == "rain") gettext("Raincloud plot") else gettext("Histogram")
    plot$plotObject <- .descriptivesSummaryScaleLayout(main, panels, mainTitle)
  }
}

# Width (px) the main plot gets in the layout, used to size its labels
.descriptivesSummaryScaleMainWidth <- function(width, panelNames) {
  if ("numbers" %in% panelNames)
    return(width * 3 / 4.75)

  if (length(panelNames) == 1) return(width / 2)

  return(width)
}

# Size of the JASP plot, depends on which panels are shown
.descriptivesSummaryScaleSize <- function(showPanels) {
  nPanels <- sum(showPanels)
  size <- switch(as.character(nPanels),
                 "0" = c(width = 720, height = 240),
                 "1" = c(width = 900, height = if (showPanels[["numbers"]]) 520 else 360),
                 "2" = c(width = 900, height = 700),
                 c(width = 1000, height = 720))

  # the numbers column needs ~330 px
  if (showPanels[["numbers"]])
    size[["width"]] <- 1120

  # scales every layout at once
  sizeFactor <- 0.9

  return(round(size * sizeFactor))
}

# Put the main plot and the panels together with patchwork
.descriptivesSummaryScaleLayout <- function(main, panels, mainTitle) {
  nPanels <- length(panels)

  if (nPanels == 0)
    return(main)

  # lettered panel titles ("A. Raincloud plot", "B. Time series", ...)
  panelTitles <- c(timeSeries = gettext("Time series"), qq = gettext("Q-Q plot"), numbers = gettext("Statistics"))
  titles <- sprintf("%s. %s", LETTERS[seq_len(nPanels + 1)], c(mainTitle, panelTitles[names(panels)]))
  main   <- .descriptivesSummaryScaleAddTitle(main, titles[1])
  panels <- Map(.descriptivesSummaryScaleAddTitle, panels, titles[-1])

  # margins around the panels so neighbouring axes are not confused; `main` has its own margins
  panels <- lapply(panels, function(pnl) pnl + ggplot2::theme(plot.margin = ggplot2::margin(t = 14, r = 10, b = 12, l = 20)))
  if (!is.null(panels[["numbers"]]))
    panels[["numbers"]] <- panels[["numbers"]] + ggplot2::theme(plot.margin = ggplot2::margin(t = 14, r = 6, b = 12, l = 12))

  # the numbers panel is a full-height column on the right (like Minitab's summary report)
  plots <- panels[names(panels) != "numbers"]

  left <- if (length(plots) == 0) {
    main
  } else if (length(plots) == 1 && is.null(panels[["numbers"]])) {
    patchwork::wrap_plots(main, plots[[1]], widths = c(1, 1))
  } else if (length(plots) == 1) {
    patchwork::wrap_plots(main, plots[[1]], ncol = 1)
  } else {
    patchwork::wrap_plots(main, patchwork::wrap_plots(plots, ncol = 2), ncol = 1)
  }

  if (is.null(panels[["numbers"]]))
    return(left)

  patchwork::wrap_plots(left, panels[["numbers"]], widths = c(3, 1.75))
}

# Panel title; on a patchwork (histogram + box) it goes on the first plot, `+` would add it to the last
.descriptivesSummaryScaleAddTitle <- function(p, title) {
  titleTheme <- ggplot2::theme(plot.title = ggplot2::element_text(size = 15, hjust = 0.5, margin = ggplot2::margin(b = 6)))

  if (inherits(p, "patchwork")) {
    p[[1]] <- p[[1]] + ggplot2::ggtitle(title) + titleTheme
    return(p)
  }

  return(p + ggplot2::ggtitle(title) + titleTheme)
}

# Median and quartiles, same quantile type as the main table
.descriptivesSummaryScaleQuartiles <- function(variable, options) {
  q <- quantile(na.omit(variable), c(0.25, 0.5, 0.75), names = FALSE,
                type = as.integer(options[["quantilesType"]]))

  return(stats::setNames(q, c("q1", "median", "q3")))
}

# Whisker ends: the most extreme values within 1.5 * IQR of the box (as geom_boxplot)
.descriptivesSummaryScaleWhiskers <- function(variable, q) {
  valid <- na.omit(variable)
  iqr   <- q[["q3"]] - q[["q1"]]
  lo    <- min(valid[valid >= q[["q1"]] - 1.5 * iqr])
  hi    <- max(valid[valid <= q[["q3"]] + 1.5 * iqr])

  return(stats::setNames(c(lo, hi), c("lo", "hi")))
}

# Width (mm) of a label: ~0.52 * fontSize per character (measured) plus 0.6 mm spacing
.descriptivesSummaryScaleLabelMM <- function(weightedNchar, fontSize) weightedNchar * fontSize * 0.52 + 0.6

# Largest font size (mm) at which no two adjacent labels of a row overlap, between minSize and maxSize.
# `hjust` decides how much of each label's width extends into the gap to its neighbour.
.descriptivesSummaryScaleFontSize <- function(vals, labelText, hjust, mainWidthPx, axisRange,
                                              minSize = 3.0, maxSize = 6) {
  n <- length(vals)
  if (n < 2) return(maxSize)

  usableWidthMM <- (mainWidthPx / 96) * 25.4 * 0.9  # 96 px/in
  mmPerUnit     <- usableWidthMM / axisRange

  fits <- sapply(seq_len(n - 1), function(i) {
    gapMM <- (vals[i + 1] - vals[i]) * mmPerUnit
    # share of each label's width that extends into this gap
    share1 <- 1 - hjust[i]
    share2 <- hjust[i + 1]
    weightedNchar <- nchar(labelText[i]) * share1 + nchar(labelText[i + 1]) * share2
    # .descriptivesSummaryScaleLabelMM() == gapMM solved for fontSize, with a 1.15 safety margin
    (gapMM - 0.6) / (weightedNchar * 0.52 * 1.15)
  })

  return(max(minSize, min(maxSize, fits)))
}

# Font size, hjust and position of the 5 labels (lower whisker, Q1, median, Q3, upper whisker), each on its true
# value. Labels are centred; a whisker label that would run past the axis edge is aligned inwards instead.
# When one row is too crowded:
# - staggered: two rows (whiskers and median on row 1, Q1 and Q3 on row 2)
# - leader lines: if even that is too tight, one row, labels moved apart with a line to their true value
# Returns list(text, size, hjust, row, pos, leader); `pos` is where each label is drawn (data units).
.descriptivesSummaryScaleLabelSizes <- function(q, w, mainWidthPx, axisRange, axisMin, axisMax,
                                                staggerBelow = 4, staggerMax = 4.5, leaderBelow = 3.5, leaderSize = 3.5) {
  vals <- c(w[["lo"]], q, w[["hi"]])
  text <- sprintf("%.2f", vals)

  usableWidthMM <- (mainWidthPx / 96) * 25.4 * 0.9
  mmPerUnit     <- usableWidthMM / axisRange

  # smallest size over the rows; minSize = 0 so a size that does not fit is not clamped up
  fit <- function(row, hjust) {
    min(sapply(split(seq_len(5), row), function(idx)
      .descriptivesSummaryScaleFontSize(vals[idx], text[idx], hjust[idx], mainWidthPx, axisRange, minSize = 0)))
  }
  # align a whisker label inwards when centred it would run past the axis edge
  clipHjust <- function(size) {
    halfWidthMM <- function(i) nchar(text[i]) * size * 0.52 / 2
    c(if (halfWidthMM(1) > (vals[1] - axisMin) * mmPerUnit) 0 else 0.5, 0.5, 0.5, 0.5,
      if (halfWidthMM(5) > (axisMax - vals[5]) * mmPerUnit) 1 else 0.5)
  }

  hjust <- rep(0.5, 5)
  row   <- rep(1, 5)
  size  <- fit(row, hjust)

  if (size < staggerBelow) {
    row  <- c(1, 2, 1, 2, 1)
    size <- min(staggerMax, fit(row, hjust))
  }

  hjust <- clipHjust(size)
  if (any(hjust != 0.5)) size <- min(size, fit(row, hjust))

  if (size >= leaderBelow || all(row == 1))
    return(list(text = text, size = rep(max(3, size), 5), hjust = hjust, row = row, pos = vals, leader = FALSE))

  # leader lines: move the labels apart with the least total displacement (isotonic regression)
  widthMM <- nchar(text) * leaderSize * 0.52 + 1.5
  c       <- c(0, cumsum((widthMM[-1] + widthMM[-5]) / 2))
  trueMM  <- (vals - axisMin) * mmPerUnit
  centre  <- stats::isoreg(seq_len(5), trueMM - c)$yf + c

  # keep the row inside the axis
  centre <- centre + max(0, widthMM[1] / 2 - centre[1])
  centre <- centre - max(0, centre[5] + widthMM[5] / 2 - axisRange * mmPerUnit)

  return(list(text = text, size = rep(leaderSize, 5), hjust = rep(0.5, 5), row = rep(1, 5),
              pos = axisMin + centre / mmPerUnit, leader = TRUE))
}

# Raincloud plot: density, box plot and optional raindrops
.descriptivesSummaryScaleRaincloud <- function(variable, options, mainWidthPx = 576) {
  valid <- na.omit(variable)
  df    <- data.frame(y = valid, x = factor(1))
  q     <- .descriptivesSummaryScaleQuartiles(variable, options)
  w     <- .descriptivesSummaryScaleWhiskers(variable, q)
  base  <- jaspGraphs::JASPcolors(options[["colorPalette"]])[1]

  showDrops <- options[["summaryScaleRaindrops"]]

  # box in the colour of the cloud: lighter fill, darker outline
  boxFill   <- grDevices::colorRampPalette(c(base, "white"))(11)[7]
  boxColour <- grDevices::colorRampPalette(c(base, "black"))(11)[5]

  # the box plot draws the outliers only when the raindrops are off, otherwise they would be drawn twice
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggrain::geom_rain(
      rain.side    = "r",
      violin.args  = list(fill = base, alpha = .6, colour = "black", linewidth = .3),
      # nudged further out than the default to leave room for the labels between box and cloud
      violin.args.pos = list(width = 0.5, position = ggplot2::position_nudge(x = 0.27)),
      boxplot.args = list(fill = boxFill, colour = boxColour, outlier.shape = if (showDrops) NA else 19, linewidth = .4),
      boxplot.args.pos = list(width = 0.05, position = ggplot2::position_nudge(x = 0.08)),
      # raindrops are drawn below instead, to colour the outliers
      point.args   = list(alpha = 0))

  if (showDrops) {
    pointsDf <- data.frame(x = factor(1), y = valid, outlier = valid < w[["lo"]] | valid > w[["hi"]])
    # same jitter as ggrain's default
    p <- p +
      ggplot2::geom_point(data = pointsDf, ggplot2::aes(x = x, y = y, colour = outlier),
                          position = ggplot2::position_jitter(width = 0.04, height = 0, seed = 42),
                          alpha = .6, size = 1.5, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`FALSE` = base, `TRUE` = "#D55E00"))
  }

  # after coord_flip() the values are still on the y scale
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(valid)

  # whisker and quartile labels in the gap between the box and the cloud
  axisRange <- diff(range(yBreaks))
  sizes <- .descriptivesSummaryScaleLabelSizes(q, w, mainWidthPx, axisRange, min(yBreaks), max(yBreaks))
  # hjust, not vjust: coord_flip() does not rotate the text itself
  rowX  <- if (sizes$leader) 1.19 else if (all(sizes$row == 1)) 1.16 else c(1.14, 1.215)[sizes$row]
  labs  <- data.frame(x = rowX, y = sizes$pos, value = c(w[["lo"]], q, w[["hi"]]), label = sizes$text,
                      size = sizes$size, hjust = sizes$hjust, face = c("plain", "plain", "bold", "plain", "plain"))
  if (sizes$leader)
    p <- p + ggplot2::geom_segment(data = labs, ggplot2::aes(x = 1.11, xend = 1.165, y = value, yend = y),
                                   colour = "grey50", linewidth = .3, inherit.aes = FALSE)
  p <- p + ggplot2::geom_label(data = labs,
                               ggplot2::aes(x = x, y = y, label = label, size = size, hjust = hjust, fontface = face),
                               fill = "white", alpha = .8, linewidth = 0, inherit.aes = FALSE) +
    ggplot2::scale_size_identity()

  # small expansions, otherwise coord_flip(xlim) leaves an empty band under the raindrops
  p <- p +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.02, 0.02))) +
    ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks), expand = ggplot2::expansion(mult = 0.02)) +
    # clip = "off" so a label at the axis edge is not cut off
    ggplot2::coord_flip(xlim = c(if (showDrops) 0.92 else 1.05, 1.54), clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    jaspGraphs::geom_rangeframe(sides = "b") +
    jaspGraphs::themeJaspRaw() +
    # left/right margins so the outer tick labels are not clipped
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(l = 14, r = 14, t = 2, b = 2))

  return(p)
}

# Histogram (Minitab style) with optional density line, box plot with labels below it on the same x axis
.descriptivesSummaryScaleHistogram <- function(variable, options, mainWidthPx = 576) {
  valid  <- variable[!is.na(variable)]
  base   <- jaspGraphs::JASPcolors(options[["colorPalette"]])[1]
  q      <- .descriptivesSummaryScaleQuartiles(variable, options)
  w      <- .descriptivesSummaryScaleWhiskers(variable, q)
  breaks <- jaspGraphs::getPrettyAxisBreaks(valid)

  h <- jaspGraphs::jaspHistogram(valid, xName = "", binWidthType = "sturges",
                                 density = options[["summaryScaleDensline"]], densityLineWidth = .8)
  h$layers[[1]]$aes_params$fill   <- base
  h$layers[[1]]$aes_params$alpha  <- .6
  h$layers[[1]]$aes_params$colour <- "black"
  h <- h +
    ggplot2::scale_x_continuous(breaks = breaks, limits = range(breaks), expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  out <- data.frame(x = valid[valid < w[["lo"]] | valid > w[["hi"]]])
  box <- data.frame(xmin = w[["lo"]], xlower = q[["q1"]], xmiddle = q[["median"]],
                    xupper = q[["q3"]], xmax = w[["hi"]], y = 0)

  boxFill   <- grDevices::colorRampPalette(c(base, "white"))(11)[7]
  boxColour <- grDevices::colorRampPalette(c(base, "black"))(11)[5]

  axisRange <- diff(range(breaks))
  sizes <- .descriptivesSummaryScaleLabelSizes(q, w, mainWidthPx, axisRange, min(breaks), max(breaks))
  labs  <- data.frame(x = sizes$pos, value = c(w[["lo"]], q, w[["hi"]]),
                      y = if (sizes$leader) 0.55 else c(0.34, 0.7)[sizes$row], label = sizes$text, size = sizes$size,
                      hjust = sizes$hjust, face = c("plain", "plain", "bold", "plain", "plain"))

  b <- ggplot2::ggplot() +
    ggplot2::geom_boxplot(data = box,
                          ggplot2::aes(y = y, xmin = xmin, xlower = xlower, xmiddle = xmiddle, xupper = xupper, xmax = xmax),
                          stat = "identity", orientation = "y", fill = boxFill, colour = boxColour,
                          linewidth = .4, width = .5) +
    ggplot2::geom_point(data = out, ggplot2::aes(x = x, y = 0), colour = boxColour, size = 1.2) +
    (if (sizes$leader) ggplot2::geom_segment(data = labs, ggplot2::aes(x = value, xend = x, y = 0.27, yend = 0.5),
                                             colour = "grey50", linewidth = .3)) +
    ggplot2::geom_text(data = labs,
                       ggplot2::aes(x = x, y = y, label = label, size = size, hjust = hjust, fontface = face),
                       vjust = 0) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(breaks = breaks, limits = range(breaks), expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(limits = c(-0.35, if (any(sizes$row == 2)) 1.15 else if (sizes$leader) 1.0 else 0.85)) +
    # so a label at the axis edge is not cut off
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    jaspGraphs::themeJaspRaw() +
    # top margin: the gap between histogram and box
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 8))

  # taller box panel when the labels take two rows
  return(patchwork::wrap_plots(h, b, ncol = 1, heights = c(4, if (any(sizes$row == 2)) 2.4 else 1.6)))
}

# Q-Q plot of the summary plots, with a 95% reference envelope (always on)
.descriptivesSummaryScaleQQ <- function(variable, options, ciLevel = 0.95) {
  .descriptivesQQNormPlot(variable, ciLevel = ciLevel, yName = gettext("Sample quantiles"),
                          fillColor = jaspGraphs::JASPcolors(options[["colorPalette"]])[1], bandOutline = FALSE)
}

# Normal Q-Q plot of the standardized values with an optional pointwise reference envelope (`ciLevel`, NULL for
# none). Same as the Distributions module (.ldFillQQPlot), on the standardized scale; Fox (2016), Section 3.1.3.
# Not jaspGraphs::plotQQnorm(): its band is centred on the quartile line, not on the y = x line it draws.
.descriptivesQQNormPlot <- function(variable, ciLevel = NULL, yName = gettext("Sample quantiles"),
                                    identicalAxes = FALSE, fillColor = "steelblue", bandOutline = TRUE) {
  z <- sort(as.numeric(scale(na.omit(variable))))
  n <- length(z)
  p <- stats::ppoints(n)

  theoretical <- stats::qnorm(p)
  df <- data.frame(theoretical = theoretical, sample = z)

  bandLayer <- NULL
  if (is.numeric(ciLevel)) {
    # standard error of the i-th order statistic, envelope centred on y = x
    slope    <- stats::coef(stats::lm(z ~ theoretical))[[2]]
    se       <- sqrt(p * (1 - p) / n) * slope / stats::dnorm(theoretical)
    crit     <- stats::qnorm((1 - ciLevel) / 2, lower.tail = FALSE)
    df$lower <- theoretical - crit * se
    df$upper <- theoretical + crit * se
    bandLayer <- ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = fillColor, alpha = .25,
                                      colour = if (bandOutline) "black" else NA)
  }

  allY    <- c(z, df$lower, df$upper)
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(if (identicalAxes) c(theoretical, allY) else theoretical)
  yBreaks <- if (identicalAxes) xBreaks else jaspGraphs::getPrettyAxisBreaks(allY)

  line <- data.frame(x = range(xBreaks))

  ggplot2::ggplot(df, ggplot2::aes(x = theoretical)) +
    bandLayer +
    ggplot2::geom_line(data = line, ggplot2::aes(x = x, y = x), colour = "darkred", linewidth = 1) +
    jaspGraphs::geom_point(ggplot2::aes(y = sample)) +
    ggplot2::scale_x_continuous(name = gettext("Theoretical quantiles"), breaks = xBreaks, limits = range(xBreaks)) +
    ggplot2::scale_y_continuous(name = yName, breaks = yBreaks, limits = range(yBreaks)) +
    jaspGraphs::geom_rangeframe() + jaspGraphs::themeJaspRaw()
}

# Time series: the values against the row number. `rowNumbers` = row.names(dataset), which keeps the original
# row numbers when a filter is active.
.descriptivesSummaryScaleTimeSeries <- function(variable, options, variableName, rowNumbers = seq_along(variable)) {
  n    <- length(variable)
  tMax <- max(rowNumbers)
  q    <- .descriptivesSummaryScaleQuartiles(variable, options)
  w    <- .descriptivesSummaryScaleWhiskers(variable, q)
  # outliers by the same 1.5 * IQR rule as the box plot
  df <- data.frame(t = rowNumbers, y = variable, outlier = variable < w[["lo"]] | variable > w[["hi"]])
  outlierColour <- "#D55E00"

  # - n <= 100: line and points at full size
  # - 100 < n <= 500: points and line shrink with n
  # - n > 500: thin line with a running median on top, only the outliers as points
  # A missing value leaves a gap in the line.
  showDots <- isTRUE(options[["summaryScaleTimeSeriesDots"]])

  p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = y))

  if (n <= 500) {
    pointSize <- if (n <= 100) 2.5 else 2.5 * sqrt(100 / n)
    # the line thins faster than the points; .12 is the width of the n > 500 line
    lineWidth <- max(.12, .4 * min(1, 100 / n))
    p <- p + ggplot2::geom_line(colour = "black", linewidth = lineWidth, na.rm = TRUE)
    if (showDots)
      p <- p +
        ggplot2::geom_point(ggplot2::aes(fill = outlier), shape = 21, colour = "black", size = pointSize,
                            stroke = .4 * pointSize / 2.5, na.rm = TRUE, show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = c(`FALSE` = "grey", `TRUE` = outlierColour))
  } else {
    # running median (robust to outliers), window ~n/20; runmed() needs an odd window and no NAs
    valid <- df[!is.na(df$y), ]
    k     <- max(3, 2 * floor(nrow(valid) / 40) + 1)
    valid$trend <- stats::runmed(valid$y, k = k, endrule = "median")

    p <- p +
      ggplot2::geom_line(colour = "grey60", linewidth = .12, na.rm = TRUE) +
      ggplot2::geom_line(data = valid, ggplot2::aes(x = t, y = trend),
                         colour = jaspGraphs::JASPcolors(options[["colorPalette"]])[1], linewidth = 1)
    if (showDots)
      p <- p + ggplot2::geom_point(data = valid[valid$outlier, ], shape = 21, fill = outlierColour, colour = "black",
                                   size = 1.5, stroke = .3)
  }

  # outliers labelled with their row number, at most `maxLabels` (the most extreme ones)
  maxLabels <- 10
  outliers  <- df[!is.na(df$outlier) & df$outlier, ]
  outliers  <- outliers[order(abs(outliers$y - q[["median"]]), decreasing = TRUE), ]
  nOutliers <- nrow(outliers)
  labelled  <- outliers[seq_len(min(nOutliers, maxLabels)), ]
  labelled$high <- labelled$y > q[["median"]]

  yBreaks <- jaspGraphs::getPrettyAxisBreaks(na.omit(variable))
  nudge   <- diff(range(yBreaks)) * .06

  for (high in c(TRUE, FALSE)) {
    side <- labelled[labelled$high == high, ]
    if (nrow(side) == 0)
      next
    p <- p + ggrepel::geom_text_repel(data = side, ggplot2::aes(x = t, y = y, label = t),
                                      nudge_y = if (high) nudge else -nudge, size = 3.5, colour = outlierColour,
                                      segment.colour = "grey50", segment.size = .3, min.segment.length = .2,
                                      box.padding = .15, point.padding = .1, seed = 1, max.overlaps = Inf)
  }

  if (nOutliers > maxLabels)
    p <- p + ggplot2::labs(subtitle = gettextf("Row numbers of the %1$i most extreme of %2$i outliers", maxLabels, nOutliers)) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, margin = ggplot2::margin(b = 4)))

  # the axis ends at the last row; drop the pretty break before it when the two labels would overlap
  tBreaks <- pretty(c(0, tMax), n = 4)
  tBreaks <- tBreaks[tBreaks < tMax]
  if (length(tBreaks) > 1 && tMax - tBreaks[length(tBreaks)] < (tBreaks[2] - tBreaks[1]) / 2)
    tBreaks <- tBreaks[-length(tBreaks)]
  tBreaks <- c(tBreaks, tMax)

  p <- p +
      ggplot2::scale_x_continuous(breaks = tBreaks, limits = c(0, tMax)) +
      # room for the nudged labels, ggplot2 drops labels outside the limits
      ggplot2::scale_y_continuous(breaks = yBreaks,
                                  limits = range(yBreaks) + c(if (any(!labelled$high)) -2 * nudge else 0,
                                                              if (any(labelled$high))  2 * nudge else 0)) +
      ggplot2::labs(x = gettext("Row number"), y = variableName) +
      jaspGraphs::geom_rangeframe() + jaspGraphs::themeJaspRaw()

  return(p)
}

# Numbers panel
.descriptivesSummaryScalePanel <- function(variable, options) {
  valid    <- na.omit(variable)
  nTotal   <- length(variable)
  nValid   <- length(valid)
  nMissing <- nTotal - nValid

  ordQuantile <- function(p) quantile(valid, p, names = FALSE, type = as.integer(options[["quantilesType"]]))

  # "-" when a statistic cannot be computed
  fmt  <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "-"
  fmtP <- function(p) if (!is.finite(p)) "-" else if (p < .001) "< .001" else sprintf("%.3f", p)

  # shapiro.test() needs 3-5000 values that are not all identical
  shapiro <- try(stats::shapiro.test(valid), silent = TRUE)
  swW <- if (isTryError(shapiro)) NaN else shapiro$statistic[[1]]
  swP <- if (isTryError(shapiro)) NaN else shapiro$p.value

  # estimate followed by its 95% CI, e.g. "43.922 [43.341, 44.503]"
  ci <- .descriptivesSummaryScaleCIs(valid)
  withCI <- function(est, ci) if (all(is.finite(ci))) sprintf("%s [%.3f, %.3f]", fmt(est), ci[1], ci[2]) else fmt(est)

  # no min and max, the main plot already shows the range
  rows <- data.frame(
    group = c(1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4),
    label = c(gettext("Total N"), gettext("Valid N"), gettext("Missing N"), gettext("Valid %"),
              gettext("Mean"), gettext("SD"),
              gettext("Skewness"), gettext("Kurtosis"), gettext("Shapiro-Wilk"), gettext("p-value"),
              gettext("25% quantile"), gettext("Median"), gettext("75% quantile")),

    value = c(nTotal, nValid, nMissing, sprintf("%.1f%%", nValid / nTotal * 100),
              withCI(mean(valid), ci[["mean"]]), withCI(sd(valid), ci[["sd"]]),
              fmt(.descriptivesSkewness(valid)), fmt(.descriptivesKurtosis(valid)), fmt(swW), fmtP(swP),
              fmt(ordQuantile(0.25)), withCI(ordQuantile(0.5), ci[["median"]]), fmt(ordQuantile(0.75))),
    stringsAsFactors = FALSE
  )

  p <- .descriptivesSummaryPanelDraw(rows, textSize = 4.5) +
    ggplot2::labs(subtitle = gettext("[ ] = 95% CI")) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, margin = ggplot2::margin(b = 4)))

  return(p)
}

# 95% CIs (as Minitab), independent of the main table's CI options. NA where they cannot be computed.
# - mean: t interval; sd: chi-square interval (n >= 2)
# - median: distribution-free from the order statistics, coverage at least 95% (n >= 6)
.descriptivesSummaryScaleCIs <- function(valid) {
  n  <- length(valid)
  na <- c(NA_real_, NA_real_)

  if (n < 2 || stats::sd(valid) == 0)
    return(list(mean = na, sd = na, median = na))

  m  <- mean(valid)
  s  <- stats::sd(valid)
  meanCI <- m + c(-1, 1) * stats::qt(.975, df = n - 1) * s / sqrt(n)
  sdCI   <- sqrt((n - 1) * s^2 / stats::qchisq(c(.975, .025), df = n - 1))

  # smallest l with P(Bin(n, .5) <= l) >= .025, so P(X(l) <= median <= X(n - l + 1)) >= .95
  l <- stats::qbinom(.025, n, .5)
  sorted   <- sort(valid)
  medianCI <- if (l >= 1) sorted[c(l, n - l + 1)] else na

  return(list(mean = meanCI, sd = sdCI, median = medianCI))
}

# "label | value" panel with lines between the groups; `rows` has the columns group, label and value
.descriptivesSummaryPanelDraw <- function(rows, textSize = 3.5) {
  rows$y <- -(seq_len(nrow(rows)) + (rows$group - 1) * 0.6)

  seps <- tapply(rows$y, rows$group, min)
  seps <- seps[-length(seps)] - 0.6

  ggplot2::ggplot() +
    ggplot2::geom_text(data = rows, ggplot2::aes(x = 0, y = y, label = label), hjust = 0, size = textSize, fontface = "bold") +
    ggplot2::geom_text(data = rows, ggplot2::aes(x = 1, y = y, label = value), hjust = 1, size = textSize) +
    ggplot2::geom_hline(yintercept = seps, linewidth = .3) +
    ggplot2::scale_x_continuous(limits = c(-0.04, 1.04)) +
    ggplot2::scale_y_continuous(limits = c(min(rows$y) - 0.7, max(rows$y) + 0.7)) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .3))
}

# Labelled empty box (placeholder, currently unused)
.descriptivesSummaryScalePlaceholder <- function(label) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = label, size = 5) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .3))
}

.geometricMean <- function(x) {
  if (any(x <= 0)) return(NA)

  return(
    exp(mean(log(x)))
  )
}

.harmonicMean <- function(x) {
  if (any(x <= 0)) return(NA)

  return(
    1 / mean(1/x)
  )
}
