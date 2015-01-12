### Make_PAR_table.R arg1 arg2 arg3
## run this script before running Model_HD_Division_Rate
## arg 1= db.name
## arg 2= output/directory
## arg 3= cruise

library(popcycle)

args <- commandArgs(TRUE)
db.name <- as.character(args[1])
out.dir <- as.character(args[2])
cruise <- as.character(args[3])

# db.name = "/Volumes/gwennm/popcycle/sqlite/popcycle.db"
# out.dir = "/Volumes/gwennm/DeepDOM/Cell_Division"
# cruise = "DeepDOM"

	sfl <- get.sfl.table(db.name)
	Par <- sfl[,c("date", "par")]
	Par$time <- strptime(Par$date, "%Y-%m-%dT%H:%M:%S", tz="GMT")
write.csv(Par[,c("par", "time")], paste0(out.dir, "/PAR_", cruise))