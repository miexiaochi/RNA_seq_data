# 加载必要的R包
library(ggplot2)
library(reshape2)

# 创建数据框，包含不同发育阶段的SOX2和NANOG表达量
expression_data <- data.frame(
  Stage = factor(rep(c("8-cell", "16-cell", "32-cell"), each = 4), levels = c("8-cell", "16-cell", "32-cell")),
  SOX2 = c(14,14,7,6,16,10,4,3,30,144,27,19),
  NANOG = c(611,277,369,214,288,239,6,183,32,477,4565,1313)
)

# 将数据转换为长格式，方便绘图
expression_data_long <- melt(expression_data, id.vars = "Stage")

# 创建箱线图
p <- ggplot(expression_data_long, aes(x = Stage, y = value, fill = variable)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16, outlier.size = 2) +
  theme_minimal() +
  labs(x = "Stage",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12), # 放大图例文本
        legend.title = element_text(size = 14), # 放大图例标题
        legend.position = "top", # 将图例放置在图像上方
        legend.box = "horizontal") + # 图例在水平方向上排列
  scale_fill_manual(values = c("SOX2" = "darkblue", "NANOG" = "coral")) +
  theme(legend.title = element_blank())
p
# 保存图像，设置图像比例为3:5
ggsave("expression_plot.pdf", plot = p, width = 10, height = 15, units = "in")
