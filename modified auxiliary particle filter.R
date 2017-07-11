library(xlsx)
library(psych)
library(MASS)
library(mnormt)
#读取数据
setwd("C:/Users/ASUS/Desktop/南方电网/协动性/数据")
y=read.xlsx("submit.xls",2,header=TRUE,stringsAsFactors=FALSE)
normalized_y=matrix(0,nrow=nrow(y),ncol=ncol(y))
#标准化
i=1
for(i in 1:ncol(y))
{
  normalized_y[,i]=(y[,i]-min(y[,i]))/(max(y[,i])-min(y[,i]))*40+60
}

#计算各指标增速，并将非指数的指标的第一个季度设为1
y_speed=y
i=1
j=1
for(i in 1:4)
{
    y_speed[1,i]=1
}

i=1
j=1
for(i in 2:nrow(y))
{
  for(j in 1:(ncol(y)-1))#最后一列是指数
  {
    y_speed[i,j]=(y[i,j]-y[i-1,j])/y[i-1,j]
  }
}
#检验多少个主成分合适，将其结果赋值给fa，其中fa.parallel的参数包括：
#自变量（去除第一列时间列），即进行主成分分析的目标；
#显示检验的项目fa，此处为“pc”（主成分）个数；
#n.iter为模拟分析迭代次数，此处设为100，可以赋值为100到1000；
#main为分析图的标题；show.legend为图例
iteration=100
legend=TRUE
main="Scree plot with parallel analysis"
fa=fa.parallel(y_speed,fa="fa",n.iter=iteration,show.legend=legend,main=main)

y_speed_mean=matrix(0,nrow=1,ncol=ncol(y))
j=1
for(j in 1:ncol(y))
{
  y_speed_mean[,j]=mean(y_speed[,j])
}
#去除均值处理
i=1
j=1
for(i in 1:nrow(y))
{
  for(j in 1:ncol(y))
  {
    y_speed[i,j]=y_speed[i,j]-y_speed_mean[,j]
  }
}


#初始化参数
x=0.1 #实际状态初始化
time=44 #44个时间点
N=2000 #粒子数
x_parameter=3#转移方程的3个参数
x_factor=1
x_estimation=matrix(0,nrow=time,ncol=1) #状态估计值
x_estimation[1,]=x #第一个状态估计值
#x_Noise=1 #系统过程噪声协方差
#x_Relation=1 #量测协方差
#z_out=matrix(0,nrow=1,ncol=time) #系统量测值矩阵
#z_out[,1]=0.1*x+0.2*x^2+0.09*x^3+rnorm(1,mean=0,sd=sqrt(x_Relation)) #第一个观测值
#x_out=matrix(0,nrow=1,ncol=time) #系统真实状态矩阵
#x_out[,1]=x #第一个真实状态

x_P=matrix(0.1,nrow=x_factor,ncol=N) #N个粒子
x_P_update=matrix(0.1,nrow=x_factor,ncol=N) #下一阶段的粒子的预测值
y_update=matrix(0.1,nrow=ncol(y),ncol=N) #通过粒子预测值得到的量测预测值
P_weight=matrix(1/N,nrow=1,ncol=N) #第二阶段权重，粒子的权重

#辅助粒子滤波变量
miu=matrix(0.1,nrow=x_factor,ncol=N)#粒子预测值的均值（或其他统计量也可以）
k=matrix(0.1, nrow = 1,ncol = N)#辅助索引，高似然程度的粒子的上标
lamda=matrix(0.1, nrow = 1,ncol = N)#辅助粒子滤波的第一阶段权重


#正则化辅助粒子滤波变量
B_P=matrix(0.1,nrow=x_parameter*x_factor,ncol=N)#参数粒子
#B_P_update=matrix(0.1,nrow=ncol(y)*5,ncol=N)
M=matrix(0.1,nrow=x_parameter*x_factor,ncol=N)#参数粒子核位置
B_P_mean=matrix(0.1,nrow=time,ncol=x_parameter*x_factor)#参数粒子均值
B_P_variance=matrix(0.1,nrow=time,ncol=x_parameter*x_factor)#参数粒子方差
a=(3*0.95-1)/(2*0.95)#收缩因子
h=sqrt(1-a)

#生成初始随机粒子
i=1
j=1
for(i in 1:x_factor)
{
  for(j in 1:N)
  {
   x_P[i,j]=x+x^2+rgamma(1, shape = 2.5,scale = 1)#在初始状态附近生成服从正态分布的粒子群
   #x_P[i,j]=x+x^2+rnorm(1, mean = 0,sd = 1)
   #x_P[i,j]=x-rgamma(1, shape = 2.5,scale = 1)
   # M[1,i]=M[1,i]
   # M[2,i]=M[2,i]
   # M[3,i]=M[3,i]
  }  
}

i=1
j=1
for(i in 1:x_parameter)
{
  for(j in 1:N)
  {
    B_P[i,j]=M[i,j]+rnorm(1,mean = 0,sd=1)
  }
}


# t=2
# for(t in 2:time)
# {
# 
#   #实际状态与观测值更新
#   x=0.5*x+2*log(abs(x))+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise))
#   #x=0.5*x+1/(1-x)+rnorm(1,mean=0,sd=sqrt(x_Noise))
#   z=0.05*x+0.05*x^2+0.05*x^3+rnorm(1,mean=0,sd=sqrt(x_Relation))
#   x_out[,t]=x
#   z_out[,t]=z
# }

t=2
for(t in 2:time)
{

  # #实际状态与观测值更新
  # x=0.5*x+25*x/(1+x^2)+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise))
  #x=0.5*x+25*x/(1+x^2)+2*log(abs(x))+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise)) 
  #x=0.5*x+2*log(abs(x))+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise))
  # #x=0.5*x+1/(1-x)+rnorm(1,mean=0,sd=sqrt(x_Noise))
  #z=0.1*x+0.2*x^2+0.09*x^3+rnorm(1,mean=0,sd=sqrt(x_Relation))
  #x_out[,t]=x
  #z_out[,t]=z
  
  i=1
  j=1
  for(i in 1:x_parameter)
  {
    for(j in 1:N)
    {
      B_P_mean[t,i]=B_P_mean[t,i]+P_weight[,j]*B_P[i,j]
    }
  }
  
  i=1
  j=1
  for(i in 1:x_parameter)
  {
    for(j in 1:N)
    {
      M[i,j]=a*B_P[i,j]+(1-a)*B_P_mean[t,i]
    }
  }
  
  i=1
  j=1
  for(i in 1:x_parameter)
  {
    for(j in 1:N)
    {
      B_P_variance[t,i]=B_P_variance[t,i]+P_weight[,j]*((B_P[i,j]-B_P_mean[t,i])^2)
    }
  }
  
  j=1
  for(j in 1:N)
  {
    miu[,j]=M[1,j]+M[2,j]*x_P[,j]+M[3,j]*x_P[,j]^2
    #miu[,j]=M[1,j]+M[2,j]*x_P[,j]
  }

  j=1
  for(j in 1:N)
  {
    #计算第一阶段权重，与似然函数p（y_t|miu_t(i)）以及各个粒子的权重成比例
    lamda[,j]=dmnorm(as.matrix(y_speed[t,]),mean=t(as.matrix(rep(miu[,j],ncol(y)))),varcov=diag(ncol(y)))*P_weight[,j]
  }

  
  #通过第一阶段权重做重要性采样，得到高似然程度的粒子的上标

  k=t(as.matrix(sample(1:N,size=N,replace = TRUE,prob = lamda)))
  
  
  
  
  i=1
  j=1
  for(i in 1:x_parameter)
  {
    for(j in 1:N)
    {
      B_P[i,j]=rnorm(1,mean=M[i,k[,j]],sd=h*sqrt(B_P_variance[t,i]))
      # B_P_update[1,i]=B_P[1,i]
      # B_P_update[2,i]=B_P[2,i]
      # B_P_update[3,i]=B_P[3,i]
    }
  }
  #计算高似然程度粒子的预测值，包括随机化过程，并计算粒子预测值得到的观测值预测值

  j=1

  for(j in 1:N)
  {
    x_P_update[,j]=B_P[1,k[,j]]+B_P[2,k[,j]]*x_P[,k[,j]]+B_P[3,k[,j]]*x_P[,k[,j]]^2+rgamma(1,shape = 2.5,scale = 1)
    #x_P_update[,j]=B_P[1,k[,j]]+B_P[2,k[,j]]*x_P[,k[,j]]-rgamma(1,shape = 2.5,scale = 1)
    #x_P_update[,j]=B_P[1,k[,j]]+B_P[2,k[,j]]*x_P[,k[,j]]+rnorm(1,mean = 0,sd = 1)
  }
  i=1
  j=1
  
  for(i in 1:ncol(y))
  {
    for(j in 1:N)
    {
      y_update[i,j]=x_P_update[,j]
    }
  }
  
  #通过高似然程度粒子预测值得到的观测值预测值与高似然粒子的下一阶段粒子期望值得到的观测值预测值的比作为各个粒子的第二阶段的权重
  j=1
    for(j in 1:N)
    {
      P_weight[,j]=dmnorm(as.matrix(y_speed[t,]),mean=t(as.matrix(y_update[,j])),varcov=diag(ncol(y)))/dmnorm(as.matrix(y_speed[t,]),mean=t(as.matrix(rep(miu[,k[,j]],ncol(y)))),varcov=diag(ncol(y)))
    }

  
  #归一化权重

  j=1
  for(j in 1:N)
  {
    P_weight[,j]=P_weight[,j]/rowSums(P_weight)
  }
  
  
  #根据第二阶段权重再次做重要性采样，得到当前阶段的粒子x_t(i)

    x_P=t(data.matrix(sample(x_P_update,size=N,replace = TRUE,prob = P_weight))) 

  #B_P[1,]=t(data.matrix(sample(B_P_update[1,],size=N,replace = TRUE,prob = P_weight))) 
  #B_P[2,]=t(data.matrix(sample(B_P_update[2,],size=N,replace = TRUE,prob = P_weight))) 
  #B_P[3,]=t(data.matrix(sample(B_P_update[3,],size=N,replace = TRUE,prob = P_weight)))
  # i=1
  # for(i in 1:N)
  # {
  #   P_weight[,i]=1/N
  # }
  #通过粒子估计状态
  x_estimation[t,]=mean(x_P) 
  
}


plot(x =seq(1:time),y=x_estimation,type = "l",col="red")

write.xlsx(x_estimation,file ="x_estimation.xls",sheetName = "sheet2",append = TRUE)
write.xlsx(B_P_mean,file ="parameter_estimation.xls",sheetName = "sheet5",append = TRUE)
lines(x_out[,1:time],col="blue")



plot(x =seq(1:time),y=x_out,col="red")
points(x_estimation[,1:time],col="blue")

sd(x_estimation)
sd(x_out)

mean(x_estimation)
mean(x_out)