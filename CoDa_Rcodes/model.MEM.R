####################################
# install MortalityForecast package
####################################

require(devtools)
devtools::install_github("mpascariu/MortalityForecast")

model_MEM = model.MEM(data = , x = 0: , y = , n = 4)
predict_model_MEM = predict(model_MEM, h = fh, x.h = 0:)
