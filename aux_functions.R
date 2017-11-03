require(plotly)
require(rbokeh)


vlnplotly<- function (object, features.plot) 
{
  data.use <- data.frame(FetchData(object = object, vars.all = features.plot), check.names = F)
  figure(data=log2(data.use))%>% ly_hexbin(nUMI,nGene,xbins = 60)
}
