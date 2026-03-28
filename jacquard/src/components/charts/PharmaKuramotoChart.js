import { useEffect, useRef, useState } from "react";

const PharmaKuramotoChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetch("/data/pharma_kuramoto_validation.json")
      .then((r) => r.json())
      .then(setData);
  }, []);

  useEffect(() => {
    if (!data || !svgRef.current) return;

    import("d3").then((d3) => {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const width = 800;
      const height = 400;
      const margin = { top: 40, right: 40, bottom: 50, left: 70 };
      const plotW = width - margin.left - margin.right;
      const plotH = height - margin.top - margin.bottom;

      const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

      // Filter drugs with both K_agg and Kuramoto_R
      const drugs = data.per_drug.filter((d) => d.K_agg != null && d.Kuramoto_R != null);

      // Also show drugs with K_agg but no R (plot at bottom)
      const drugsKOnly = data.per_drug.filter((d) => d.K_agg != null && d.Kuramoto_R == null);

      const allWithK = data.per_drug.filter((d) => d.K_agg != null);

      const xScale = d3.scaleLog()
        .domain([50, 500000])
        .range([0, plotW]);

      const yScale = d3.scaleLinear()
        .domain([0.5, 1.0])
        .range([plotH, 0]);

      // Color by class
      const classes = [...new Set(allWithK.map((d) => d.class))];
      const colorScale = d3.scaleOrdinal()
        .domain(classes)
        .range(["#58E6D9", "#B63E96", "#f59e0b", "#22c55e", "#8b5cf6", "#ef4444"]);

      // Grid
      g.append("g").selectAll("line")
        .data(yScale.ticks(5))
        .join("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d))
        .attr("stroke", "#2a2a2a").attr("stroke-dasharray", "2,4");

      // Axes
      g.append("g")
        .attr("transform", `translate(0,${plotH})`)
        .call(d3.axisBottom(xScale).ticks(5, ".0e"))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      g.append("g")
        .call(d3.axisLeft(yScale).ticks(5))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      // Therapeutic threshold line at R=0.7
      g.append("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", yScale(0.7)).attr("y2", yScale(0.7))
        .attr("stroke", "#f59e0b")
        .attr("stroke-width", 1.5)
        .attr("stroke-dasharray", "8,4");

      g.append("text")
        .attr("x", plotW - 5).attr("y", yScale(0.7) - 6)
        .attr("text-anchor", "end")
        .attr("fill", "#f59e0b")
        .attr("font-size", "10px")
        .text("Therapeutic Threshold R=0.7");

      // Drug points with Kuramoto R
      drugs.forEach((d) => {
        g.append("circle")
          .attr("cx", xScale(d.K_agg))
          .attr("cy", yScale(d.Kuramoto_R))
          .attr("r", 8)
          .attr("fill", colorScale(d.class))
          .attr("stroke", "#1b1b1b")
          .attr("stroke-width", 2)
          .attr("opacity", 0.9);

        g.append("text")
          .attr("x", xScale(d.K_agg) + 12)
          .attr("y", yScale(d.Kuramoto_R) + 4)
          .attr("fill", "#e5e5e5")
          .attr("font-size", "10px")
          .attr("font-weight", "500")
          .text(d.name);
      });

      // Drugs with K_agg only (shown as triangles at bottom)
      drugsKOnly.forEach((d) => {
        const tx = xScale(d.K_agg);
        const ty = plotH - 10;
        g.append("polygon")
          .attr("points", `${tx},${ty - 6} ${tx - 5},${ty + 4} ${tx + 5},${ty + 4}`)
          .attr("fill", colorScale(d.class))
          .attr("opacity", 0.6);

        g.append("text")
          .attr("x", tx)
          .attr("y", ty + 16)
          .attr("text-anchor", "middle")
          .attr("fill", "#999")
          .attr("font-size", "8px")
          .text(d.name);
      });

      // Legend
      const legend = g.append("g").attr("transform", `translate(10, 5)`);
      classes.forEach((cls, i) => {
        legend.append("circle")
          .attr("cx", 6).attr("cy", i * 18 + 6)
          .attr("r", 5).attr("fill", colorScale(cls));
        legend.append("text")
          .attr("x", 16).attr("y", i * 18 + 10)
          .attr("fill", "#ccc").attr("font-size", "10px").text(cls);
      });

      // Title
      g.append("text").attr("x", plotW / 2).attr("y", -18)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "14px")
        .attr("font-weight", "600").text("Pharmacological Kuramoto Validation");

      // Axis labels
      g.append("text").attr("x", plotW / 2).attr("y", plotH + 42)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px")
        .text("Aggregate Coupling K_agg (log scale)");

      g.append("text").attr("transform", "rotate(-90)").attr("x", -plotH / 2).attr("y", -50)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px")
        .text("Kuramoto Order Parameter R");
    });
  }, [data]);

  return (
    <svg
      ref={svgRef}
      viewBox="0 0 800 400"
      className="w-full h-auto"
      style={{ background: "#1b1b1b", borderRadius: "8px" }}
    />
  );
};

export default PharmaKuramotoChart;
