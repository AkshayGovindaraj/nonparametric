
trial <- matrix(c(20,30,30,20), ncol=2)
colnames(trial) <- c('Democrat', 'Republican')
rownames(trial) <- c('Male', 'Female')
party.table <- as.table(trial)

barplot_graph <- barplot(party.table,beside=T,legend=T)
Chisquare_table<- chisq.test(party.table,correct = T) 
attributes(Chisquare_table)
Chisquare_table$p.value 

trial <- matrix(c(25,6,8,20), ncol=2)
colnames(trial) <- c('Pass', 'Fail')
rownames(trial) <- c('Attended', 'Skipped')
grade.table <- as.table(trial)

barplot_graph <- barplot(grade.table,beside=T,legend=T)
Fisher_test<- fisher.test(grade.table) # p-value = 7.428e-05 < 0.05. Hence we reject the null hypotheses.

trial <- matrix(c(189,10845,104,10933), ncol=2)
colnames(trial) <- c('Placebo ', 'Aspirin')
rownames(trial) <- c('Myocardial Infarction', 'No MI')
aspirin.table <- as.table(trial)

barplot_graph <- barplot(aspirin.table,beside=T,legend=T)
mcnemar.test(aspirin.table)
