#初始化参数
x=0.1 #实际状态初始化
x_Noise=1 #系统过程噪声协方差
x_Relation=1 #量测协方差
time=75 #75个时间点
N=1000 #粒子数

z_out=matrix(0,nrow=1,ncol=time) #实际观测值矩阵
z_out[,1]=x^2/20+rnorm(1,mean=0,sd=sqrt(x_Relation)) #第一个观测值

x_out=matrix(0,nrow=1,ncol=time) #真实状态矩阵
x_out[,1]=x #第一个真实状态
x_estimation=matrix(0,nrow=1,ncol=time) #状态估计值
x_estimation[,1]=x #第一个状态估计值

x_P_update=matrix(0,nrow=1,ncol=N) #下一阶段粒子的预测值
z_update=matrix(0,nrow=1,ncol=N) #通过粒子预测值得到的观测预测值
P_weight=matrix(1/N,nrow=1,ncol=N) #第二阶段权重，粒子的权重

#辅助粒子滤波变量
lamda=matrix(0, nrow = 1,ncol = N)#辅助粒子滤波的第一阶段权重
k=matrix(0, nrow = 1,ncol = N)#辅助索引，按照似然程度为概率值抽得的粒子的上标
miu=matrix(0,nrow=1,ncol=N)#下一阶段粒子预测值的均值（或其他统计量）

#初始化分布（正态）
Variance=2 #初始方差
x_P=matrix(0,nrow=1,ncol=N) #N个粒子

#生成初始随机粒子
i=1
for(i in 1:N)
{
  x_P[,i]=x+rnorm(1,mean=0,sd=sqrt(Variance)) #在初始状态附近生成服从正态分布的粒子群
}#上述循环中的正态分布即初始的后验概率分布

t=2
for(t in 2:time)
{

  #实际状态与观测值更新
  x=0.5*x+25*x/(1+x^2)+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise)) 
  z=x^2/20+rnorm(1,mean=0,sd=sqrt(x_Relation))
  x_out[,t]=x
  z_out[,t]=z 

  i=1
  for(i in 1:N)
  {
     #计算下一阶段粒子的期望,也即没有随机化的x_P_update，也可以计算其它统计量，不过均是与p{x_t|x_t-1(i)}相关的统计量
     miu[,i]=0.5*x_P[,i]+25*x_P[,i]/(1+x_P[,i]^2)+8*cos(1.2*(t-1))

     #计算第一阶段权重，与似然函数p（y_t|miu_t(i)）以及各个粒子的权重成比例
     lamda[,i]=dnorm(z,mean=miu[,i]^2/20,sd=x_Relation)*P_weight[,i]

  }
  #通过第一阶段权重做重要性采样，得到高似然程度的粒子的上标，也就是通过统计量miu估计的观测值与实际值越接近，其被抽中的概率越大
  k=t(data.matrix(sample(1:N,size=N,replace = TRUE,prob = lamda))) 
  
  #计算高似然程度粒子的预测值，即加入噪声，并计算粒子预测值得到的观测值预测值
  i=1
  for(i in 1:N)
  {
    x_P_update[,i]=0.5*x_P[,k[,i]]+25*x_P[,k[,i]]/(1+x_P[,k[,i]]^2)+8*cos(1.2*(t-1))+rnorm(1,mean=0,sd=sqrt(x_Noise))
    z_update[,i]=x_P_update[,i]^2/20
   }
  
  #通过高似然程度粒子预测值得到的观测值预测值与高似然粒子的粒子期望值得到的观测值预测值的比作为各个粒子的第二阶段的权重，注意z_update已经是通过
  #x_P[,k[,i]]得到的预测值了，所以不需要写作z_update[,k[,i]]
  i=1
  for(i in 1:N)
  {
    P_weight[,i]=dnorm(z,mean=z_update[,i],sd=x_Relation)/dnorm(z,mean=miu[,k[,i]]^2/20,sd=x_Relation)
  }

  #归一化权重
  i=1
  for(i in 1:N)
  {
    P_weight[,i]=P_weight[,i]/rowSums(P_weight)
  }

  #根据第二阶段权重再次做重要性采样，得到当前阶段的粒子x_t(i)
  x_P=t(data.matrix(sample(x_P_update,size=N,replace = TRUE,prob = P_weight))) 
  
  #通过粒子估计状态
  x_estimation[,t]=mean(x_P[1,]) 
}



plot(x =seq(1:time),y=x_estimation,type="l",col="red")
lines(x_out[,1:time],col="blue")


