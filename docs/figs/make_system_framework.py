# -*- coding: utf-8 -*-
"""
生成"双平台协同纯方位多目标跟踪系统"框架流程图。

输出: docs/figs/system_framework.{png,svg,pdf}

用法:
    python docs/figs/make_system_framework.py

字体: 标题/框名 用黑体(SimHei), 说明文字 用宋体(SimSun), 公式用 STIX。
脚本内置"文字越界"检查: 任何文字超出其方框会打印 [OVERFLOW] 提示。
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.font_manager import FontProperties

plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["svg.fonttype"] = "path"      # SVG 文字转路径，跨机显示一致
plt.rcParams["pdf.fonttype"] = 42          # PDF 内嵌 TrueType 字体

F_TITLE = FontProperties(family="SimHei", size=13.0)
F_SUB = FontProperties(family="SimSun", size=7.8)
F_NAME = FontProperties(family="SimHei", size=9.0)
F_DETAIL = FontProperties(family="SimSun", size=7.6)
F_CONT = FontProperties(family="SimHei", size=8.8)
F_NOTE = FontProperties(family="SimSun", size=7.0)
F_SMALL = FontProperties(family="SimHei", size=8.6)

W, H = 82.0, 98.0
COL = "#44546a"
C_ENV = dict(fc="#eef3f8", ec="#5a7d9a")
C_MOD = dict(fc="#eaf0f7", ec="#40658f")
C_ASS = dict(fc="#e9f4ec", ec="#4f8a5b")
C_FUS = dict(fc="#f2edf8", ec="#72559b")
C_MAN = dict(fc="#fdf1e3", ec="#b07a3a")
C_ECM = dict(fc="#fbeaea", ec="#a05252")
C_EVA = dict(fc="#f2f2f2", ec="#777777")

SCALE = 0.088                            # 1 逻辑单位对应的英寸数 (决定成图物理尺寸)
fig = plt.figure(figsize=(W * SCALE, H * SCALE))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.set_aspect("equal")
ax.axis("off")

BOXES = []


def box(x0, y0, w, h, name, detail="", style=C_MOD, lw=1.15, name_fp=None,
        detail_fp=None, pad=1.35, gap=2.0, r=0.95, tag=""):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=style["ec"],
                                facecolor=style["fc"], zorder=3))
    cx = x0 + w / 2.0
    ny = y0 + h - pad
    t1 = ax.text(cx, ny, name, ha="center", va="top", zorder=4,
                 fontproperties=name_fp or F_NAME, color="#12202e")
    t2 = None
    if detail:
        t2 = ax.text(cx, ny - gap, detail, ha="center", va="top", zorder=4,
                     fontproperties=detail_fp or F_DETAIL, color="#2b3a48",
                     linespacing=1.45)
    BOXES.append(dict(tag=tag or name, x=(x0, x0 + w), y=(y0, y0 + h),
                      texts=[t for t in (t1, t2) if t is not None]))


def container(x0, y0, w, h, title, align="left"):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h,
                                boxstyle="round,pad=0,rounding_size=1.3",
                                linewidth=0.95, edgecolor="#9bb0c4",
                                facecolor="none", linestyle=(0, (5, 3)), zorder=1))
    if align == "center":
        ax.text(x0 + w / 2.0, y0 + h - 0.5, title, ha="center", va="top", zorder=2,
                fontproperties=F_CONT, color="#3d5a76")
    else:
        ax.text(x0 + 0.9, y0 + h - 0.5, title, ha="left", va="top", zorder=2,
                fontproperties=F_CONT, color="#3d5a76")


def _fancy(p0, p1):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=9.5,
                                 linewidth=1.1, color=COL, shrinkA=0, shrinkB=0,
                                 zorder=5))


def arrow(p0, p1):
    _fancy(p0, p1)


def poly_arrow(points):
    for a, b in zip(points[:-2], points[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=COL, lw=1.1, zorder=5,
                solid_capstyle="round")
    _fancy(points[-2], points[-1])


# ---------------------------------------------------------------- 标题
ax.text(W / 2, 95.2, "双平台协同纯方位多目标跟踪系统框架", ha="center", va="center",
        fontproperties=F_TITLE, color="#0f1b28")
ax.text(W / 2, 91.9, "IMM-UKF 滤波　·　JPDA 关联　·　三角交会初始化　·　协方差交叉融合　·　平台机动优化",
        ha="center", va="center", fontproperties=F_SUB, color="#4a5a6a")
ax.plot([3, W - 3], [90.2, 90.2], color="#c9d4df", lw=0.9, zorder=1)

# ---------------------------------------------------------------- 输入层
box(2.5, 84.4, 12.5, 5.0, "场景定义", "scenario_defs.m",
    style=C_ENV, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=6.9))
box(19.0, 84.4, 26.0, 5.0, "目标真值轨迹", "generate_truth.m", style=C_ENV,
    name_fp=F_SMALL, detail_fp=FontProperties(family="SimSun", size=6.9))
arrow((15.2, 86.9), (18.8, 86.9))

box(7.0, 73.2, 38.0, 8.4, "双站方位量测",
    "generate_measurements.m\n"
    r"$z_{i,k}=\arctan2(\Delta y_{i,k},\Delta x_{i,k})+\eta_{i,k}$"
    "\n"
    r"$P_d$ 检测　$\lambda_c$ 杂波　频谱特征（可选）",
    style=C_ENV, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=7.0), gap=1.9)
box(50.0, 73.2, 22.0, 8.4, "观测平台状态",
    r"$p_{w1,k},\,p_{w2,k}$" "\n" "ECM 决定下一步位置\n" "6 m/s，仅航向可变",
    style=C_ENV, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=7.0), gap=1.9)

arrow((26.0, 84.2), (26.0, 81.8))
arrow((49.8, 77.4), (45.2, 77.4))

# ---------------------------------------------------------------- IMM 预测
container(5.5, 60.2, 63.5, 10.6, "① IMM 预测（每条航迹独立执行）")
box(7.5, 61.0, 18.0, 7.4, "输入混合",
    r"$c_j=\sum_i\Pi_{ij}\mu_i$" "\n" r"$\mu_{i|j}=\Pi_{ij}\mu_i/c_j$",
    style=C_MOD, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size= 6.9), gap=1.55)
box(28.0, 61.0, 18.0, 7.4, "三模型时间更新",
    "CV / CT / CA" "\n" r"$P\leftarrow P+Q$",
    style=C_MOD, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=6.9), gap=2.0)
box(48.5, 61.0, 20.0, 7.4, "组合量测预测",
    r"$z_{comb},\ S_{comb}$" "\n" "含模型间离散项",
    style=C_MOD, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=6.9), gap=2.0)

arrow((25.7, 64.7), (27.8, 64.7))
arrow((46.2, 64.7), (48.3, 64.7))
arrow((26.0, 73.0), (26.0, 70.9))

# ---------------------------------------------------------------- 关联 + 更新
container(5.5, 37.0, 63.5, 19.0, "② 双平台顺序 JPDA 关联 + IMM-UKF 更新", align="center")

box(8.0, 48.2, 27.0, 5.6, "平台 1：JPDA 关联",
    r"波门 → 联合事件枚举 → $\beta^{(1)}_t(j),\ \beta^{(1)}_0$",
    style=C_ASS, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=7.2), gap=1.8)
box(40.0, 48.2, 27.0, 5.6, "平台 2：JPDA 关联",
    r"波门 → 联合事件枚举 → $\beta^{(2)}_t(j),\ \beta^{(2)}_0$",
    style=C_ASS, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=7.2), gap=1.8)

box(8.0, 38.0, 27.0, 7.6, "逐模型状态更新",
    r"$x_j\leftarrow x_j+K_j\bar\nu_j$" "\n"
    r"$P_j\leftarrow P_j-(1-\beta_0)K_jS_jK_j^{T}$",
    style=C_MOD, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size=7.2), gap=2.0)
box(40.0, 38.0, 27.0, 7.6, "逐模型状态更新",
    r"$\mu_j\propto c_j\,N(\bar\nu_j;0,S_j)$" "\n"
    r"$x=\sum_j\mu_jx_j$",
    style=C_MOD, name_fp=F_SMALL,
    detail_fp=FontProperties(family="SimSun", size= 7.2), gap=1.55)

arrow((21.5, 48.0), (21.5, 45.8))
arrow((53.5, 48.0), (53.5, 45.8))
arrow((16.5, 60.0), (16.5, 54.0))
arrow((62.0, 60.0), (62.0, 54.0))

# ---------------------------------------------------------------- CI 融合
box(10.0, 27.2, 52.0, 7.6, "③ 双平台信息融合（逐模型 CI）",
    r"$P^{-1}=P_1^{-1}+P_2^{-1}-P_{pred}^{-1}$" "\n"
    r"$\mu\propto\sqrt{\mu_1\mu_2}$ 归一化后按 $\mu$ 组合状态",
    style=C_FUS, name_fp=FontProperties(family="SimHei", size=9.2),
    detail_fp=FontProperties(family="SimSun", size= 7.4), gap=1.7)

arrow((21.5, 37.8), (31.5, 35.0))
arrow((53.5, 37.8), (43.5, 35.0))

# ---------------------------------------------------------------- 航迹管理
box(6.0, 16.4, 60.0, 8.6, "④ 航迹管理　manage_tracks.m",
    r"删除：$n_{miss}\geq3$ 或越出监视区　　确认：$n_{hit}\geq3$" "\n"
    r"出生：未关联量测三角交会 → 区域 / 重复 / 频谱校验 → IMM 初始化",
    style=C_MAN, name_fp=FontProperties(family="SimHei", size=9.2),
    detail_fp=FontProperties(family="SimSun", size=7.4), gap=2.6)

arrow((36.0, 27.0), (36.0, 25.2))

# ---------------------------------------------------------------- 记录评价
box(18.0, 8.0, 38.0, 7.0, "记录与评价",
    "tracks_history → 轨迹 / 方位 / RMSE" "\n"
    "GOSPA（loc / missed / false）",
    style=C_EVA, name_fp=FontProperties(family="SimHei", size=9.0),
    detail_fp=FontProperties(family="SimSun", size=7.4), gap=2.2)

arrow((34.0, 16.2), (34.0, 15.2))

# ---------------------------------------------------------------- ECM 回路
box(69.5, 43.0, 11.5, 8.4, "ECM 平台机动",
    "贪心搜索\n7 个候选航向\n" r"$\pm30^\circ$",
    style=C_ECM, name_fp=FontProperties(family="SimHei", size=8.4),
    detail_fp=FontProperties(family="SimSun", size=6.9), gap=1.9)

poly_arrow([(66.2, 20.7), (75.2, 20.7), (75.2, 42.8)])
poly_arrow([(75.2, 51.6), (75.2, 77.4), (72.2, 77.4)])
ax.text(77.0, 33.0, r"$k\leftarrow k+1$", rotation=90, ha="center", va="center",
        fontproperties=FontProperties(family="SimSun", size=8.0), color="#4a5a6a")

# ---------------------------------------------------------------- 脚注
note_lines = [
    "注：实线箭头为数据流，虚线框为模块分组；框内的 .m 为对应实现文件。",
    r"参数：$P_d=0.9$，$\lambda_c=2$ 帧$^{-1}$站$^{-1}$，$\sigma_\beta=1^\circ$，"
    r"$\sigma_a^{CV}=0.5$ m/s$^2$，模型集 $\{$CV, CT, CA$\}$，$\Pi_{ii}=0.94$。",
    r"说明：平台 2 在未含平台 1 更新的预测副本上独立关联；协方差更新含散布项 $K\sigma_\nu^2K^{T}$；"
    r"ECM 目标 $\max\sum_t\log(1+H_tP_tH_t^{T}/R)$。",
    r"状态 $x=[p_x,p_y,v_x,v_y,\omega,a_x,a_y]^{T}$：$\omega$ 由 CT 模型估计，$(a_x,a_y)$ 由 CA 模型估计。",
]
ax.text(2.5, 1.0, "\n".join(note_lines), ha="left", va="bottom",
        fontproperties=FontProperties(family="SimSun", size=6.4),
        color="#4a5a6a", linespacing=1.55)

# ---------------------------------------------------------------- 越界自检
fig.canvas.draw()
ren = fig.canvas.get_renderer()
bad = 0
for b in BOXES:
    x0, x1 = b["x"]
    y0, y1 = b["y"]
    for t in b["texts"]:
        bb = t.get_window_extent(ren).transformed(ax.transData.inverted())
        if bb.x0 < x0 - 0.15 or bb.x1 > x1 + 0.15 or bb.y0 < y0 - 0.15:
            bad += 1
            print(f"[OVERFLOW] {b['tag']}: text x=({bb.x0:.1f},{bb.x1:.1f}) "
                  f"box x=({x0:.1f},{x1:.1f})  text y0={bb.y0:.1f} box y0={y0:.1f}")
print("overflow items:", bad)

here = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(here, f"system_framework.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white")
print("saved png / svg / pdf")
