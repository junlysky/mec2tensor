# 震源机制及矩张量参数转换的计算小工具
## 功能
### 1. mec2ten： 
#### 假设震源模型为双力偶模型，由震源机制节面（任意1个节面的strike,dip,rake）解计算地震矩张量的全部参数。
#### 使用方法：编译mec2ten.f90后，执行二进制程序
#### 输入文件名词为：mec.d
#### mec.d的数据格式：lon lat strike dip rake
#### 执行程序后，输出文件为：mecinfo.dat  

### 2. pt2ten：  
#### 假设震源模型为双力偶模型，由PT轴信息计算震源机制的节面解及其对应矩张量参数。
### 3. tensor2mec：  
#### 假设震源模型为双力偶模型，由地震矩张量计算节面解及PT轴参数。 
### 4. 
#### 计算单个震源机制参数，并绘制震源球
