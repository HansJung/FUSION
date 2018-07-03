suppressMessages(library('foreach'))
suppressMessages(library('Matrix'))
suppressMessages(library('splines'))
suppressMessages(library('SuperLearner'))
suppressMessages(library('xgboost'))
suppressMessages(library('gam'))
suppressMessages(library('gbm'))
suppressMessages(library('glmnet'))
suppressMessages(library('randomForest'))

# assumption
# treatment: singleton, binary
# outcome: singleton, binary
estimateCausalEffect <- function (file, prediction, predictionDiscrete, numFolds, numIntervals, task) {
    dataFrame = fileToDataFrame(file)
    algorithm = selectAlgorithm(dataFrame, task, predictionDiscrete)
    result = estimate(dataFrame, algorithm, prediction, numFolds, numIntervals, task)
    
    return(result)
}

selectAlgorithm <- function (dataFrame, task, predictionDiscrete) {
    if (isBinary(task[['treatment']]))
        return(predictionDiscrete)
    else
        return('linear')
}

estimate <- function (dataFrame, algorithm, prediction, numFolds, numIntervals, task) {
    if (algorithm == 'linear')
        return(estimateLinear(dataFrame, prediction, numFolds, numIntervals, task))
    else if (algorithm == 'AIPW')
        return(estiamteAIPW(dataFrame, prediction, numFolds, task))
    else if (algorithm == 'TMLE')
        return(estimateTMLE(dataFrame, prediction, numFolds, task))
}

estimateLinear <- function (dataFrame, prediction, numFolds, numIntervals, task) {
    treatment = task[['treatment']]

    X = getDataFrameColumn(dataFrame, treatment, TRUE)

    numSamples = nrow(dataFrame)

    if (isDiscrete(treatment))
        intervals = sort(getUniqueValues(dataFrame, getLabels(treatment)))
    else
        intervals = getIntervals(X, numIntervals)

    predCausal = predictLinearByInterval(dataFrame, prediction, numFolds, task, intervals)
    predObs = predictLinearByInterval(dataFrame, prediction, numFolds, task, intervals, TRUE)

    dataCausal = getSummaryLinear(predCausal, numSamples)
    dataObs = getSummaryLinear(predObs, numSamples)

    result = list()
    result[['interval']] = intervals
    result[['causal']] = dataCausal
    result[['observational']] = dataObs

    return(result)
}

predictLinearByInterval <- function (dataFrame, prediction, numFolds, task, intervals, isObservational = FALSE) {
    treatment = task[['treatment']]
    outcome = task[['outcome']]
    variables = task[['variables']]

    X = getDataFrameColumn(dataFrame, treatment, TRUE)
    Y = getDataFrameColumn(dataFrame, outcome, TRUE)
    Z = getDataFrameColumn(dataFrame, variables)
    XZ = dataFrame[ , -which(names(dataFrame) %in% getLabels(outcome))]

    numSamples = nrow(dataFrame)
    doX = matrix(rep(intervals, each=numSamples), nrow=numSamples)

    if (isObservational) {
        Yx = c()

        linearModel = trainModel(Y, data.frame(X), prediction, numFolds, outcome)

        for (i in 1:length(intervals)) {
            Xi = cbind(doX[, i])
            colnames(Xi) = names(X)
            predict = predictModel(linearModel, Xi)
            Yx = cbind(Yx, predict)
        }

        return(Yx)
    } else {
        hat_Yx = c()

        linearModel = trainModel(Y, XZ, prediction, numFolds, outcome)

        for (i in 1:length(intervals)) {
            XiZ = cbind(doX[, i], Z)
            colnames(XiZ) = names(XZ)
            predict = predictModel(linearModel, XiZ)
            hat_Yx = cbind(hat_Yx, predict)
        }

        return(hat_Yx)
    }
}

estimateIPW <- function (dataFrame, prediction, numFolds, task) {
    treatment = task[['treatment']]
    outcome = task[['outcome']]
    variables = task[['variables']]

    # avoid dividing by 0
    epsilon = 1e-8

    X = getDataFrameColumn(dataFrame, treatment, TRUE)
    Y = getDataFrameColumn(dataFrame, outcome, TRUE)
    Z = getDataFrameColumn(dataFrame, variables)

    propScoreModel = trainModel(X, Z, prediction, numFolds, outcome)
    propScorePrediction = predictModel(propScoreModel, Z)

    pX1GivenZ = propScorePrediction + epsilon
    pX0GivenZ = (1 - pX1GivenZ) + epsilon

    Y_X0_reweighted = (Y / pX0GivenZ) * (1 - X)
    Y_X1_reweighted = (Y / pX1GivenZ) * X

    reweighted = Y_X0_reweighted + Y_X1_reweighted

    reweightDataset = cbind(dataFrame)
    reweightDataset[, getLabels(outcome)] = reweighted
    reweightDatasetString = paste('', apply(reweightDataset, 1, paste, collapse = ' '))

    result = list()
    result[['0']] = pX0GivenZ
    result[['1']] = pX1GivenZ
    result[['reweighted']] = reweightDatasetString
    result[['predict']] = propScorePrediction

    return(result)
}

estiamteAIPW <- function (dataFrame, prediction, numFolds, task) {
    treatment = task[['treatment']]
    outcome = task[['outcome']]
    variables = task[['variables']]

    X = getDataFrameColumn(dataFrame, treatment, TRUE)
    Y = getDataFrameColumn(dataFrame, outcome, TRUE)
    XZ = dataFrame[ , -which(names(dataFrame) %in% getLabels(outcome))]
    X0Z = getInterventionalDistribution(dataFrame, treatment, 0)
    X1Z = getInterventionalDistribution(dataFrame, treatment, 1)
    
    linearModel = trainLinearModel(Y, XZ, numFolds, outcome)
    hat_Yx0 = predictModel(linearModel, X0Z)
    hat_Yx1 = predictModel(linearModel, X1Z)

    propScoreResult = estimateIPW(dataFrame, prediction, numFolds, task)
    pX0GivenZ = propScoreResult[['0']]
    pX1GivenZ = propScoreResult[['1']]

    Y_X0 = (Y / pX0GivenZ) * (1 - X) - ((1 - X) - pX0GivenZ) / (pX0GivenZ) * mean(hat_Yx0)
    Y_X1 = (Y / pX1GivenZ) * X - ((X) - pX1GivenZ) / (pX1GivenZ) * mean(hat_Yx1)

    # observational
    XY = dataFrame[ , -which(names(dataFrame) %in% getLabels(variables))]
    Y_X0_obs = dataFrame[dataFrame[getLabels(treatment)] == 0, 1]
    Y_X1_obs = dataFrame[dataFrame[getLabels(treatment)] == 1, 1]

    datasetCausal = cbind(Y_X0, Y_X1)
    datasetObs = cbind(Y_X0_obs, Y_X1_obs)

    result = list()
    result[['interval']] = c(0,1)
    result[['causal']] = getSummaryBinaryTreatment(datasetCausal)
    result[['observational']] = getSummaryBinaryTreatment(datasetObs)
    result[['reweighted']] = propScoreResult[['reweighted']]

    return(result)
}

estimateTMLE <- function (dataFrame, prediction, numFolds, task) {
    treatment = task[['treatment']]
    outcome = task[['outcome']]
    variables = task[['variables']]

    X = getDataFrameColumn(dataFrame, treatment, TRUE)
    Y = getDataFrameColumn(dataFrame, outcome, TRUE)
    XZ = dataFrame[ , -which(names(dataFrame) %in% getLabels(outcome))]
    X0Z = getInterventionalDistribution(dataFrame, treatment, 0)
    X1Z = getInterventionalDistribution(dataFrame, treatment, 1)

    linearModel = trainLinearModel(Y, XZ, numFolds, outcome)
    hat_Yx0 = predictModel(linearModel, X0Z)
    hat_Yx1 = predictModel(linearModel, X1Z)
    hat_Y_xz = predictModel(linearModel, XZ)

    propScoreResult = estimateIPW(dataFrame, prediction, numFolds, task)
    pX0GivenZ = propScoreResult[['0']]
    pX1GivenZ = propScoreResult[['1']]

    H1 = X / pX1GivenZ
    H0 = (1 - X) / pX0GivenZ
    H = H1 + H0

    TMLE_update = lm(Y ~ -1 +offset(hat_Y_xz) + H)
    eps = TMLE_update$coef

    Y_TMLE_X0 =  hat_Yx0 + eps * H0
    Y_TMLE_X1 =  hat_Yx1 + eps * H1

    # observational
    XY = dataFrame[ , -which(names(dataFrame) %in% getLabels(variables))]
    Y_X0_obs = dataFrame[dataFrame[getLabels(treatment)] == 0, 1]
    Y_X1_obs = dataFrame[dataFrame[getLabels(treatment)] == 1, 1]

    datasetCausal = cbind(Y_TMLE_X0, Y_TMLE_X1)
    datasetObs = cbind(Y_X0_obs, Y_X1_obs)

    result = list()
    result[['interval']] = c(0,1)
    result[['causal']] = getSummaryBinaryTreatment(datasetCausal)
    result[['observational']] = getSummaryBinaryTreatment(datasetObs)
    result[['reweighted']] = propScoreResult[['reweighted']]

    return(result)
}

trainModel <- function (dependent, regressors, prediction, numFolds, outcome) {
    predictions = toSLPrediction(prediction)
    propScoreModel = SuperLearner(Y = dependent, X = regressors, family=getFamilyFunction(outcome), cvControl=list(V=numFolds), SL.library=predictions)
    return(propScoreModel)
}

predictModel <- function (model, regressors) {
    predictObject = predict(model, regressors, onlySL = T)
    return(predictObject$pred)
}

trainLinearModel <- function (dependent, regressors, numFolds, outcome) {
    linearPredictionAlgorithm = c('SL.lm')
   linearModel = SuperLearner(Y = dependent, X = regressors, family=getFamilyFunction(outcome), cvControl = list(V=numFolds), SL.library=linearPredictionAlgorithm)
    return(linearModel)
}

getSummaryLinear <- function (predictDataPoints, numSamples) {
    zValue = 1.96

    statistics = getStatistics(predictDataPoints)

    upperDataPoints = statistics[['mean']] + zValue * statistics[['sd']] / sqrt(numSamples)
    lowerDataPoints = statistics[['mean']] - zValue * statistics[['sd']] / sqrt(numSamples)

    result = list()
    result[['data']] = list()
    result[['data']][['predict']] = statistics[['mean']]
    result[['data']][['predict_upper']] = upperDataPoints
    result[['data']][['predict_lower']] = lowerDataPoints
    result[['statistics']] = statistics

    return(result)
}

getSummaryBinaryTreatment <- function (dataPoints) {
    result = list()
    result[['data']] = list()
    result[['data']][['predict']] = dataPoints
    result[['statistics']] = getStatistics(dataPoints)

    return(result)
}

getStatistics <- function (dataPoints) {
    mean = apply(dataPoints, 2, mean)
    sd = apply(dataPoints, 2, sd)

    getQ1 <- function (dataPoints) {
        summary(dataPoints)[2]
    }

    getQ3 <- function (dataPoints) {
        summary(dataPoints)[5]
    }

    q1 = apply(dataPoints, 2, getQ1)
    q3 = apply(dataPoints, 2, getQ3)

    result = list()
    result[['mean']] = mean
    result[['sd']] = sd
    result[['q1']] = q1
    result[['q3']] = q3

    return(result)
}

getInterventionalDistribution <- function (dataFrame, treatment, xValue) {
    numRows = nrow(dataFrame)

    dataFrame[, getLabels(treatment)] = as.vector(rep(xValue, numRows))

    return(dataFrame)
}

fileToDataFrame <- function (file) {
    con <- textConnection(file)
    dataFrame <- read.table(con, header=TRUE, comment.char="")
    close(con)

    return(dataFrame)
}

getDataFrameColumn <- function (dataFrame, variables, isSingleVariable = FALSE) {
    result = subset(dataFrame, select=getLabels(variables))

    if (isSingleVariable)
        return(c(t(result)))
    else
        return(result)
}

getFamilyFunction <- function (outcome) {
    if (isBinary(outcome))
        return(binomial())
    else
        return(gaussian())
}

toSLPrediction <- function (predictions) {
    SLPredictions <- c()

    for (pred in predictions) {
        if (pred == 0)
            SLPredictions <- c(SLPredictions, 'SL.xgboost')
        else if (pred == 1)
            SLPredictions <- c(SLPredictions, 'SL.randomForest')
        else if (pred == 2)
            SLPredictions <- c(SLPredictions, 'SL.cforest')
        else if (pred == 3)
            SLPredictions <- c(SLPredictions, 'SL.gbm')
        else if (pred == 4)
            SLPredictions <- c(SLPredictions, 'SL.gam')
        else if (pred == 5)
            SLPredictions <- c(SLPredictions, 'SL.glm')
        else if (pred == 6)
            SLPredictions <- c(SLPredictions, 'SL.logreg')
    }

    return(SLPredictions)
}

getIntervals <- function (dataPoints, numIntervals) {
    return(seq(min(dataPoints), max(dataPoints), length.out = numIntervals))
}

getUniqueValues <- function (dataFrame, labels) {
    return(unique(dataFrame[, labels]))
}

isDiscrete <- function (variables) {
    if (0 %in% getTypes(variables))
        return(TRUE)
    return(FALSE)
}

isBinary <- function (variables) {
    if (variables[4][1][['isBinary']])
        return(TRUE)
    return(FALSE)
}

getLabels <- function (variables) {
    return(variables[1][1][['label']])
}

getTypes <- function (variables) {
    return(variables[3][1][['type']])
}

# getPlotImage <- function (dataPoints) {
#     png(filename="temp.png", width=500, height=500)
#     plot(density(dataPoints))
#     dev.off()

#     plotBinary <- paste(readBin("temp.png", what="raw", n=1e7), collapse="")
#     return(plotBinary)
# }