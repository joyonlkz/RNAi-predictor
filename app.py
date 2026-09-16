import streamlit as st
import pandas as pd
import plotly.express as px
from Bio.Seq import Seq

# ---------- 序列处理函数 ----------
def clean_sequence(seq: str) -> str:
    seq = seq.strip().upper().replace("\n", "").replace("\r", "").replace(" ", "")
    seq = seq.replace("T", "U")
    return "".join(ch for ch in seq if ch in "AUGC")

def gc_content(seq: str) -> float:
    if not seq: return 0.0
    return round(sum(1 for b in seq if b in "GC") / len(seq) * 100, 1)

def get_antisense(sense: str) -> str:
    return str(Seq(sense).reverse_complement_rna())

def has_internal_repeat(seq: str, min_rep: int = 4) -> bool:
    return any(len(set(seq[i:i+min_rep])) == 1 for i in range(len(seq)-min_rep+1))

def count_au(seq: str) -> int:
    return sum(1 for b in seq if b in "AU")

# ---------- Reynolds评分 ----------
def score_reynolds(sense: str, antisense: str) -> dict:
    score, details = 0, []
    gc = gc_content(sense)
    if 30 <= gc <= 52: score += 1; details.append("GC合格")
    if len(sense) >= 19 and sense[18] == "A": score += 1; details.append("19位A")
    if count_au(sense[14:19]) >= 3: score += 1; details.append("15-19位A/U充足")
    if len(sense) >= 3 and sense[2] == "A": score += 1; details.append("3位A")
    if len(sense) >= 10 and sense[9] == "U": score += 1; details.append("10位U")
    if len(sense) >= 13 and sense[12] != "G": score += 1; details.append("13位非G")
    if not has_internal_repeat(sense): score += 1; details.append("无内部重复")
    if antisense and antisense[0] in "GC": score += 1; details.append("反义5'端G/C")
    return {"score": score, "details": ", ".join(details)}

# ---------- Streamlit界面 ----------
st.set_page_config(page_title="RNAi Predictor", layout="wide")
st.title("🧬 RNAi siRNA 预测工具")
st.caption("仅供研究使用，基于Reynolds规则进行候选筛选，不保证实验效果。")

mrna_input = st.text_area("输入mRNA序列（纯序列或FASTA格式）", height=200,
                          placeholder=">GeneName\nATGGTGAGCAAGGGC...")

col1, col2, col3 = st.columns(3)
with col1: sirna_len = st.selectbox("siRNA长度", [19, 21, 23], index=1)
with col2: top_n = st.slider("显示Top N", 10, 50, 20)
with col3:
    gc_min = st.number_input("GC最低%", 0, 100, 30)
    gc_max = st.number_input("GC最高%", 0, 100, 52)

if st.button("🚀 开始预测", type="primary"):
    if not mrna_input.strip():
        st.warning("请输入序列"); st.stop()
    raw = mrna_input.strip()
    seq = clean_sequence("".join(raw.split("\n")[1:]) if raw.startswith(">") else raw)
    if len(seq) < sirna_len:
        st.error(f"序列长度不足{sirna_len}nt"); st.stop()
    if len(seq) > 5000:
        st.error("序列长度不能超过5000nt"); st.stop()

    with st.spinner("正在生成候选..."):
        candidates = [{"pos": i+1, "sense": seq[i:i+sirna_len]}
                      for i in range(len(seq)-sirna_len+1)]
    results = []
    for c in candidates:
        gc = gc_content(c["sense"])
        if not (gc_min <= gc <= gc_max): continue
        anti = get_antisense(c["sense"])
        sc = score_reynolds(c["sense"], anti)
        results.append({"位置": c["pos"], "正义链": c["sense"], "反义链": anti,
                        "GC%": gc, "Reynolds得分": sc["score"], "评分详情": sc["details"]})
    if not results:
        st.warning("没有符合条件的候选，请放宽GC范围"); st.stop()

    df = pd.DataFrame(results).sort_values("Reynolds得分", ascending=False).head(top_n)
    df.insert(0, "排名", range(1, len(df)+1))
    st.success(f"共生成 {len(results)} 个候选，展示Top {len(df)}")
    st.dataframe(df, use_container_width=True, hide_index=True)

    fig = px.scatter(df, x="位置", y="Reynolds得分", color="Reynolds得分",
                     color_continuous_scale=["red","yellow","green"], range_color=[0,8],
                     hover_data=["正义链","GC%"], title="候选siRNA分布")
    st.plotly_chart(fig, use_container_width=True)

    csv = df.to_csv(index=False).encode("utf-8-sig")
    st.download_button("⬇️ 下载CSV", data=csv, file_name="rnai_candidates.csv", mime="text/csv")

with st.expander("📋 评分标准说明（Reynolds规则）"):
    st.markdown("GC 30-52% | 19位A | 15-19位≥3个A/U | 3位A | 10位U | 13位非G | 无内部重复 | 反义5'端G/C。**总分≥6判定为有效候选。**")