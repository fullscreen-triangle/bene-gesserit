import { useEffect, useRef, useState } from "react";

const TripleEquivalenceChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetch("/data/triple_equivalence_empirical.json")
      .then((r) => r.json())
      .then(setData);
  }, []);

  useEffect(() => {
    if (!data || !svgRef.current) return;

    import("d3").then((d3) => {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const width = 600;
      const height = 400;
      const margin = { top: 40, right: 30, bottom: 50, left: 70 };
      const plotW = width - margin.left - margin.right;
      const plotH = height - margin.top - margin.bottom;

      const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

      const elements = data.per_element;
      const modalities = data.modalities;

      const yScale = d3.scaleBand()
        .domain(elements.map((d) => `${d.symbol} (Z=${d.Z})`))
        .range([0, plotH])
        .padding(0.15);

      const xScale = d3.scaleBand()
        .domain(modalities)
        .range([0, plotW])
        .padding(0.15);

      // Grid
      g.append("g").selectAll("line")
        .data(elements)
        .join("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", (d) => yScale(`${d.symbol} (Z=${d.Z})`) + yScale.bandwidth())
        .attr("y2", (d) => yScale(`${d.symbol} (Z=${d.Z})`) + yScale.bandwidth())
        .attr("stroke", "#2a2a2a");

      // Cells
      elements.forEach((el) => {
        modalities.forEach((mod) => {
          const yKey = `${el.symbol} (Z=${el.Z})`;
          g.append("rect")
            .attr("x", xScale(mod))
            .attr("y", yScale(yKey))
            .attr("width", xScale.bandwidth())
            .attr("height", yScale.bandwidth())
            .attr("fill", el.all_agree ? "#22c55e" : "#ef4444")
            .attr("opacity", el.all_agree ? 0.7 : 0.7)
            .attr("rx", 4);

          g.append("text")
            .attr("x", xScale(mod) + xScale.bandwidth() / 2)
            .attr("y", yScale(yKey) + yScale.bandwidth() / 2)
            .attr("dy", "0.35em")
            .attr("text-anchor", "middle")
            .attr("fill", "#fff")
            .attr("font-size", "11px")
            .attr("font-weight", "600")
            .text(el.all_agree ? "\u2713" : "\u2717");
        });
      });

      // Y axis (element labels)
      g.append("g")
        .call(d3.axisLeft(yScale).tickSize(0))
        .call((el) => el.select(".domain").remove())
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "11px"));

      // X axis (modality labels)
      g.append("g")
        .attr("transform", `translate(0,${plotH})`)
        .call(d3.axisBottom(xScale).tickSize(0))
        .call((el) => el.select(".domain").remove())
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "11px"));

      // Title
      g.append("text").attr("x", plotW / 2).attr("y", -18)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "14px")
        .attr("font-weight", "600").text("Triple Equivalence: Four-Modality Agreement");

      // Summary
      g.append("text").attr("x", plotW / 2).attr("y", plotH + 40)
        .attr("text-anchor", "middle").attr("fill", "#22c55e").attr("font-size", "12px")
        .attr("font-weight", "600")
        .text(`${data.total_agreements}/${data.total_checks} checks passed (${(data.overall_agreement_rate * 100).toFixed(0)}%)`);
    });
  }, [data]);

  return (
    <svg
      ref={svgRef}
      viewBox="0 0 600 400"
      className="w-full h-auto"
      style={{ background: "#1b1b1b", borderRadius: "8px" }}
    />
  );
};

export default TripleEquivalenceChart;
