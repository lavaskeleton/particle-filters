#初始化参数
x=0.1 #实际状态初始化
x_Noise=1 #转移方程噪声协方差
x_Relation=1 #量测方程噪声协方差
time=75 #75个时间点
N=1000 #粒子数

z_out=matrix(0,nrow=1,ncol=time) #实际量测值矩阵
z_out[,1]=x^2/20+rnorm(1,mean=0,sd=sqrt(x_Relation)) #第一个观测值

x_out=matrix(0,nrow=1,ncol=time) #实际状态矩阵
x_out[,1]=x #第一个真实状态
x_estimation=matrix(0,nrow=1,ncol=time) #状态估计值
x_estimation[,1]=x #第一个状态估计值

x_P_update=matrix(0,nrow=1,ncol=N) #下一阶段粒子预测值，服从预测分布p{x_t|y_1:t-1}
z_update=matrix(0,nrow=1,ncol=N) #通过粒子预测值得到的量测预测值,服从分布p{y_t|x_t(i)*}
P_weight=matrix(0,nrow=1,ncol=N) #权重


#初始化分布（正态）
Variance=2 #初始方差
x_P=matrix(0,nrow=1,ncol=N) #N个粒子

#生成初始随机粒子
i=1
for(i in 1:N)
{
  x_P[,i]=x+rnorm(1,mean=0,sd=sqrt(Variance)) #在初始状态附近生成服从正态分布的粒子群
}#上述循环中的正态分布即初始的后验概率分布p(x_0)

t=2
for(t in 2:time)
{
  #实际状态与观测值更新
  x=0.5*x+25*x/(1+x^2)+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise)) 
  z=x^2/20+rnorm(1,mean=0,sd=sqrt(x_Relation))
  x_out[,t]=x#x在现实中是未知的，所以下列程序中不应该出现x或x_out
  z_out[,t]=z#可用信息 
  
  i=1
  for(i in 1:N)
  {
    #预测阶段
    ##通过状态转移方程得到下一阶段的粒子，也就是状态的预测值x_t(i)*，服从预测分布p{x_t|y_1:t-1}
    x_P_update[,i]=0.5*x_P[,i]+25*x_P[,1]/(1+x_P[,i]^2)+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise))
    ##通过新的粒子群得到观测值的预测值y_t(i)*，服从分布p{y_t|x_t(i)*}
    z_update[,i]=x_P_update[,i]^2/20
    
    #更新阶段
    ##计算粒子的未归一化的权重（使得观测值的预测值与实际值值越接近的粒子权重越大）
    P_weight[,i]=dnorm(z,mean=z_update[,i],sd=x_Relation)
  }
  ##归一化权重
  i=1
  for(i in 1:N)
  {
    P_weight[,i]=P_weight[,i]/rowSums(P_weight)
  }
  
  ##根据归一化的权重抽样，得到x_t(i)
  x_P=t(data.matrix(sample(x_P_update,size=N,replace = TRUE,prob = P_weight))) 
  
  x_estimation[,t]=mean(x_P[1,]) #通过粒子估计真实状态x_t（由于重抽样，所有粒子的权重变为1/N）
}

#作图对比准确率，一般粒子数越多越准确？
plot(x =seq(1:time),y=x_estimation,type="l",col="red")
lines(x_out[,1:time],col="blue")


