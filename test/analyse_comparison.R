library(abcrf)

npods = 200
# entry data
neutral = read.table("/home/croux/Programmes/msms/ABCsweep/test/neutral/rep1/outputABC_sumStats.txt", h=T)
neutral = rbind(neutral, read.table("/home/croux/Programmes/msms/ABCsweep/test/neutral/rep2/outputABC_sumStats.txt", h=T))
neutral = na.omit(neutral)
neutral_references = neutral[-((nrow(neutral)-npods):nrow(neutral)),]
neutral_pods = neutral[(nrow(neutral)-npods):nrow(neutral),]

selection = read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep1/outputABC_sumStats.txt", h=T)
selection = rbind(selection, read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep2/outputABC_sumStats.txt", h=T))
selection = na.omit(selection)
selection_references = selection[1:nrow(neutral_references),]
selection_pods = selection[(nrow(selection)-npods):nrow(selection),]

modIndex = as.factor(c(rep("neutral", nrow(neutral_references)), rep("selection", nrow(selection_references))))
mod = abcrf(modIndex, rbind(neutral_references, selection_references), ntree = 2000, paral=T)
print(mod)

predict_neutral = predict(mod, neutral_pods, ntree=2000, paral=T)
predict_selection = predict(mod, selection_pods, ntree=2000, paral=T)

