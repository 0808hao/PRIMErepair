# PRIMErepair 使用说明

这是第 131–133 步方法的首个 R 包版本 0.1.0。它接收锁定预测，进行候选选择与有限审计池认证；不训练模型，不读取 eICU 原始大表，也不自动判定原始特征是否存在时间泄露。

## 安装和首次运行

在 R 里运行（将源码包放到当前工作目录）：

```r
install.packages("PRIMErepair.tar.gz", repos = NULL, type = "source")
library(PRIMErepair)
ex <- prime_example()
fit <- prime_audit(ex$pool, prime_oracle(ex$y), epsilon = 0.002)
print(fit)
```

R 版本要求至少 4.3；除 R 标准包外没有额外依赖，也不需要编译器。Linux、Windows 和 macOS 都使用源码安装方式；本次实际检查的平台记录在完整压缩包的 reports 中，不能把未测试平台写成已验证。

## 接入自己的锁定预测

`predictions` 是样本 × 候选的数值概率矩阵。必须有唯一列名；每个候选与参考预测必须对应同一批样本、同一结局。`cost` 是输入成本，必须提前规定。参考候选成本须为最高（可并列）。`groups` 是行数相同的 0/1 矩阵，列是需要约束的亚组；可重叠，但不能为空。

```r
# 可执行示例；真实数据时用自己的锁定矩阵替换这里的 p。
p <- ex$pool$public$P
cost <- setNames(ex$pool$candidates$cost, ex$pool$candidates$candidate)
pool <- prime_pool(p, cost, reference = "reference",
  groups = ex$pool$public$G, ids = ex$pool$public$ids)
fit <- prime_audit(pool, prime_oracle(ex$y),
  epsilon = 0.002, budgets = c(0.1, 0.25, 0.5, 1), batch = 64)
```

有名字的成本按候选名匹配；无名字的成本按预测矩阵列顺序匹配。概率矩阵和亚组矩阵如果有行名，必须与 ids 完全同序。无名字的结局按样本行顺序；有名字的回调结果必须与本次请求的样本 ID 同序。

不要把新分析数据塞进现有 pool 后修改字段；内容校验会拒绝已被改动的对象。应重新创建 pool，并把它视为新的审计设计。

## 三种不同输出

1. `prime_audit(..., method="decision_focused")`：确定性支持界。预测、结局、成本和池定义满足约定时，界限覆盖所有未查询二分类标签的可能补全。
2. `prime_audit(..., method="decision_cell_residual_wor_EB", alpha=0.05)`：概率模式。允许预设错误概率，层内必须均匀无放回。包只改变选哪个公开分区，不按分区内个人分数挑标签。
3. `prime_lure(...)`：固定预算点估计。输出明确写 `NONE_POINT_ESTIMATE_ONLY`，不能写成认证。

两种认证模式都只保证这个固定审计池和锁定候选库中的目标，不保证未来人群或临床安全。0.002 是此前测试用的工程阈值，不是已经临床论证的不劣界值。输入成本与审计标签成本分开计量。

## 读结果

`fit$summary` 每行对应一个预算；概率模式还区分 epsilon 和 alpha。

- `selected_candidate`：当前已证实可行的候选。
- `selected_cost`：其声明输入成本。
- `certified_minimum_cost`：是否已排除所有更便宜候选。FALSE 时不能把返回的可行候选称为“已认证最低成本”。
- `labels_used`：该配置沿当前查询路径达到状态所用标签。
- `fit$queries`：实际买过的所有唯一标签；`nrow(fit$queries)` 才是整个共享多配置运行的总付费数。
- `fit$bounds`：候选 × 亚组界限。正的 Brier 差表示候选比参考差。
- `proxy_completed_mean`：仅为查询启发式，不是随机抽样无偏估计，也不是置信区间。

停止条件只要求最低成本，不要求同成本模型全部排出名次。因此不同策略返回不同但同成本且合格的模型，可以同时正确。

`budgets` 先换算为 `max(1, floor(n * fraction))` 的整数标签上限。小于一个标签的比例会变成一个标签。预算节点、批量大小、随机种子和配置顺序均属于协议；比较方法时保持一致。修改预算节点可能改变批次边界和停止标签数。

## 真实标签服务

```r
# 真正使用时，实现自己的回调。这里使用模拟标签演示契约。
reveal <- function(indices) {
  ids <- pool$public$ids[indices]
  setNames(ex$y[indices], ids)
}
fit <- prime_audit(pool, reveal)
```

包只请求尚未付费的行，每个已查询单位的标签可同时服务所有候选。若服务失败，会抛出 `prime_query_error`，其中保留 `completed_queries` 和 `pending_indices`。已经在外部收费但未返回的请求无法由本包判断；应先核对外部账目。本版不提供跨进程请求队列或持久化续跑接口。

`prime_oracle(y)` 是有状态的，每次独立审计需要重新创建。重复使用同一 oracle 会拒绝重复收费。各策略比较必须分别创建 oracle，不能让某个策略免费使用另一个策略的答案。

## 保存表格

```r
out <- file.path(tempdir(), paste0("PRIMErepair_", format(Sys.time(), "%Y%m%d_%H%M%S")))
prime_export(fit, out)
```

输出包括：

- `01_methodology_tables`：协议、候选、亚组。
- `02_result_tables`：决策、界限、实际查询标签与诊断。
- `audit_result.rds`：完整结果对象。
- `file_manifest.csv`：各文件校验值。

已存在的目录不会被覆盖。查询表含实际输入的样本 ID 和结局；不要把自己的真实查询表误当作随包示例发布。压缩包自带的可执行数据全部是模拟数据，历史 eICU/Sleep 文件仅保留汇总统计。

## 并行

单次自适应查询路径按顺序执行。独立池、独立随机种子或独立模拟重复可在外层用 R 的 `parallel::mclapply(..., mc.cores=20)` 并行；Windows 可用 PSOCK。不要在并行任务之间共享有状态的 oracle。每个任务明确自己的种子，并在任务内部创建 oracle。本包恢复调用者的 RNG 状态。

## 当前证据边界

eICU 主分析支持较低的确定性认证标签成本。Sleep 的确定性策略与 shared-width 打平，概率策略还可更省；Sleep 保留双特征族的结论本身也是有效决策。首版不承诺每个数据集都领先，不把多个医院划分称为多个独立队列，不把方法学认证写成临床安全证明。

完整包内的服务器复核脚本只重放已有锁定预测和标签，不重读 lab、不重训模型。它是部署环境的复现检查，不是另做一轮调参。
