View(data)

prog_start_end=do.call(rbind, by(data[!is.na(data$progression),], INDICES = data$id[!is.na(data$progression)], 
         function(x){
          if(tail(x$progression,1)==1){
            tail(x$time, 2)
          }else{
            c(tail(x$time,1), Inf)
          }
         }
))

npmle = interval::icfit(Surv(prog_start_end[,1],prog_start_end[,2],type="interval2")~1)
