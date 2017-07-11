library(xlsx)
library(psych)

setwd("C:/Users/ASUS/Desktop/南方电网/协动性/数据")
y=read.xlsx("submit.xls",1,header=TRUE,stringsAsFactors=FALSE)
normalized_y=matrix(0,nrow=nrow(y),ncol=ncol(y))
i=1
for(i in 1:ncol(y))
{
  normalized_y[,i]=(y[,i]-min(y[,i]))/(max(y[,i])-min(y[,i]))*40+60
}
#检验多少个主成分合适，将其结果赋值给fa，其中fa.parallel的参数包括：
#自变量（去除第一列时间列），即进行主成分分析的目标；
#显示检验的项目fa，此处为“pc”（主成分）个数；
#n.iter为模拟分析迭代次数，此处设为100，可以赋值为100到1000；
#main为分析图的标题；show.legend为图例
iteration=100
legend=TRUE
main="Scree plot with parallel analysis"
fa=fa.parallel(y,fa="fa",n.iter=iteration,show.legend=legend,main=main)







#初始化参数
x=0.1 #实际状态初始化
x_Noise=1 #系统过程噪声协方差
x_Relation=1 #量测协方差
time=44 #75个时间点
N=2000 #粒子数


z_out=matrix(0,nrow=1,ncol=time) #系统量测值矩阵
z_out[,1]=0.1*x+0.2*x^2+0.09*x^3+rnorm(1,mean=0,sd=sqrt(x_Relation)) #第一个观测值

x_out=matrix(0,nrow=1,ncol=time) #系统真实状态矩阵
x_out[,1]=x #第一个真实状态
x_estimation=matrix(0,nrow=1,ncol=time) #状态估计值
x_estimation[,1]=x #第一个状态估计值

x_P=matrix(0,nrow=1,ncol=N) #N个粒子
x_P_update=matrix(0,nrow=1,ncol=N) #下一阶段的粒子的预测值
z_update=matrix(0,nrow=1,ncol=N) #通过粒子预测值得到的量测预测值
P_weight=matrix(1/N,nrow=1,ncol=N) #第二阶段权重，粒子的权重

#辅助粒子滤波变量
miu=matrix(0.1,nrow=1,ncol=N)#粒子预测值的均值（或其他统计量也可以）
k=matrix(0.1, nrow = 1,ncol = N)#辅助索引，高似然程度的粒子的上标
lamda=matrix(0.1, nrow = 1,ncol = N)#辅助粒子滤波的第一阶段权重


#正则化辅助粒子滤波变量
B_P=matrix(0.1,nrow=3,ncol=N)
B_P_update=matrix(0.1,nrow=3,ncol=N)
M=matrix(0.01,nrow=3,ncol=N)
B_P_mean=matrix(0.1,nrow=3,ncol=time)
B_P_variance=matrix(0.1,nrow=3,ncol=time)
a=(3*0.95-1)/(2*0.95)
h=sqrt(1-a)

#初始化分布（正态）
Variance=1 #初始方差

#生成初始随机粒子
i=1
for(i in 1:N)
{
  x_P[,i]=x+rnorm(1,mean=0,sd=sqrt(Variance)) #在初始状态附近生成服从正态分布的粒子群
  # M[1,i]=M[1,i]
  # M[2,i]=M[2,i]
  # M[3,i]=M[3,i]
  B_P[1,i]=M[1,i]+rnorm(1,mean=0,sd=sqrt(Variance))
  B_P[2,i]=M[2,i]+rnorm(1,mean=0,sd=sqrt(Variance))
  B_P[3,i]=M[3,i]+rnorm(1,mean=0,sd=sqrt(Variance))
  
}#上述循环中的正态分布即初始的后验概率分布


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
  x=0.5*x+2*log(abs(x))+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise))
  # #x=0.5*x+1/(1-x)+rnorm(1,mean=0,sd=sqrt(x_Noise))
  z=0.1*x+0.2*x^2+0.09*x^3+rnorm(1,mean=0,sd=sqrt(x_Relation))
  x_out[,t]=x
  z_out[,t]=z
  
  i=1
  for(i in 1:N)
  {
    B_P_mean[1,t]=B_P_mean[1,t]+P_weight[,i]*B_P[1,i]
    B_P_mean[2,t]=B_P_mean[2,t]+P_weight[,i]*B_P[2,i]
    B_P_mean[3,t]=B_P_mean[3,t]+P_weight[,i]*B_P[3,i]
    # P_weight[,i]=1/N
  }
  i=1
  for (i in 1:N)
  {
    B_P_variance[1,t]=B_P_variance[1,t]+P_weight[,i]*((B_P[1,i]-B_P_mean[1,t])^2)
    B_P_variance[2,t]=B_P_variance[2,t]+P_weight[,i]*((B_P[2,i]-B_P_mean[2,t])^2)
    B_P_variance[3,t]=B_P_variance[3,t]+P_weight[,i]*((B_P[3,i]-B_P_mean[3,t])^2)
  }
  
  i=1
  for(i in 1:N)
  {
    #计算下一阶段参数粒子的核位置
    M[1,i]=a*B_P[1,i]+(1-a)*B_P_mean[1,t]
    M[2,i]=a*B_P[2,i]+(1-a)*B_P_mean[2,t] 
    M[3,i]=a*B_P[3,i]+(1-a)*B_P_mean[3,t] 
  }
  i=1
  j=1
  for(i in 1:N)
  {
    miu[,i]=0.5*x_P[,i]+2*log(abs(x_P[,i]))+8*cos(1.2*(t-1))
  }
  
  i=1
  for(i in 1:N)
  {
    #计算第一阶段权重，与似然函数p（y_t|miu_t(i)）以及各个粒子的权重成比例
    lamda[,i]=dnorm(z,mean=M[1,i]*miu[,i]+M[2,i]*miu[,i]^2+M[3,i]*miu[,i]^3,sd=x_Relation)*P_weight[,i]
  }
  
  #通过第一阶段权重做重要性采样，得到高似然程度的粒子的上标
  k=t(data.matrix(sample(1:N,size=N,replace = TRUE,prob = lamda))) 
  
  
  i=1
  for(i in 1:N)
  {
    B_P[1,i]=rnorm(1,mean=M[1,k[,i]],sd=h*sqrt(B_P_variance[1,t]))
    B_P[2,i]=rnorm(1,mean=M[2,k[,i]],sd=h*sqrt(B_P_variance[2,t]))
    B_P[3,i]=rnorm(1,mean=M[3,k[,i]],sd=h*sqrt(B_P_variance[3,t]))
    
    B_P_update[1,i]=B_P[1,i]
    B_P_update[2,i]=B_P[2,i]
    B_P_update[3,i]=B_P[3,i]
  }
  #计算高似然程度粒子的预测值，包括随机化过程，并计算粒子预测值得到的观测值预测值
  i=1
  for(i in 1:N)
  {
    x_P_update[,i]=0.5*x_P[,k[,i]]+2*log(abs(x_P[,k[,i]]))+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise))
    z_update[,i]=B_P[1,i]*x_P_update[,i]+B_P[2,i]*x_P_update[,i]^2+B_P[3,i]*x_P_update[,i]^3
  }
  
  
  #通过高似然程度粒子预测值得到的观测值预测值与高似然粒子的下一阶段粒子期望值得到的观测值预测值的比作为各个粒子的第二阶段的权重
  i=1
  for(i in 1:N)
  {
    P_weight[,i]=dnorm(z,mean=z_update[,i],sd=x_Relation)/dnorm(z,mean=M[1,k[,i]]*miu[,k[,i]]+M[2,k[,i]]*miu[,k[,i]]^2+M[3,k[,i]]*miu[,k[,i]]^3,sd=x_Relation)
  }
  
  #归一化权重
  i=1
  for(i in 1:N)
  {
    P_weight[,i]=P_weight[,i]/rowSums(P_weight)
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
  x_estimation[,t]=mean(x_P) 
  
}


plot(x =seq(1:time),y=x_out,type = "l",col="red")
lines(x_estimation[,1:time],col="blue")


plot(x =seq(1:time),y=x_out,col="red")
points(x_estimation[,1:time],col="blue")

sd(x_estimation)
sd(x_out)

mean(x_estimation)
mean(x_out)