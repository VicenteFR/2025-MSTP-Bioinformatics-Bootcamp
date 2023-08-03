# Library
library(fmsb)

#TJ plot

data <- t(as.data.frame(matrix(c(5,9,6,10,7,3,8,7,6))))
colnames(data) <- c("basketball","R", "python","annoying enthusiasm" , "opinion of T-Cells", "dancing skills","willingness to dance","league of legends","lobster catching")

#Hannah plot

data <- t(as.data.frame(matrix(c(10,7,8,10,9,8,4,7,4))))
colnames(data) <- c("writing grants and manuscripts","R", "python","figuring out \nwhat you did wrong \nbased on your plot" , "tumor evolution", "machine learning","eating out","crashing rafts","showing up late to meetings")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data <- rbind(rep(10,9) , rep(0,9) , data)
data<-as.data.frame(data)
# Check your data, it has to look like this!

# Custom the radarChart !
radarchart( data  , axistype=1 , 
            
            #custom polygon
            pcol=rgb(0.2,0.5,0.5,0.9) , 
            pfcol=rgb(0.2,0.5,0.5,0.5)  ,
            plwd=2 ,  caxislabels = c(2, 4, 6, 8,10) ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
            
            #custom labels
            vlcex=1.1
)

