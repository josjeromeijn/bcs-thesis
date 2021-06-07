# fit equation size_t = Linf-(Linf-Lb)*exp(-r_B*age)

input2 <- read.table("LOCATION\\lm.txt", header=TRUE,  sep = "\t")  # define correct path where file is located

id_list <- input2[[1]] #create list with idnumbers

vBgrowthrate_list <- data.frame(matrix(ncol = 4, nrow = 0)) #create outcome list to store outcomes
colnames(vBgrowthrate_list) <- c("")
names(vBgrowthrate_list) <- c("Estimate","Std. Error","t value","Pr(>|t|)")

for (val in id_list) { #iterate over idnumbers
  index <- match(c(val),id_list) #find the index of idnumber

  file_name = paste('id',val,'.txt',sep="") #create filename of idnumber
  input <- read.table(paste("LOCATION",file_name,sep=""), header=TRUE,  sep = "\t") #load in file of idnumber
  Lb = input['size_t'][1,1]; #find Lb of idnumber  
  Linf = tail(input['size_t'],n=1)[1,1] #find the Linf of idnumber
  vBgrowthrate <- nls(size_t~Linf-(Linf-Lb)*exp(-b*age),data=input, start=c(b=0.000001),control = nls.control(maxiter = 100))
  vBgrowthrate_list <- rbind(vBgrowthrate_list,summary(vBgrowthrate)$parameters) #add parameters to outcome list
  }

rownames(vBgrowthrate_list) <- id_list #change rownames to idnumbers
write.table(vBgrowthrate_list, "LOCATION\\growthrates.txt",sep="\t")

 