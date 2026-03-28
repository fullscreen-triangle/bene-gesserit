import { useEffect, useRef, useState } from "react";

const ComputingStackChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetch("/data/full_computing_stack_validation.json")
      .then((r) => r.json())
      .then(setData);
  }, []);

  useEffect(() => {
    if (!data || !svgRef.current) return;

    import("d3").then((d3) => {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const width = 800;
      const height = 300;
      const margin = { top: 40, right: 180, bottom: 30, left: 200 };
      const plotW = width - margin.left - margin.right;
      const plotH = height - margin.top - margin.bottom;

      const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

      const layers = data.layers;

      // Build metric labels
      const metricLabels = layers.map((layer) => {
        switch (layer.layer) {
          case 0: return `Rectification: ${layer.rectification_ratio.toFixed(1)}`;
          case 1: return `On/Off: ${layer.on_off_ratio}`;
          case 2: return `Gates: ${layer.average_agreement_pct}%`;
          case 3: return `Fidelity: ${(layer.average_fidelity * 100).toFixed(1)}%`;
          case 4: return `Speedup: ${layer.speedup_factor.toLocaleString()}x`;
          case 5: return `Reliability: ${(layer.fibonacci_reliability * 100).toFixed(0)}%`;
          default: return "";
        }
      });

      // Normalize to 0-100 for display
      const barValues = layers.map((layer) => {
        switch (layer.layer) {
          case 0: return Math.min(layer.rectification_ratio, 100);
          case 1: return Math.min(layer.on_off_ratio, 100);
          case 2: return layer.average_agreement_pct;
          case 3: return layer.average_fidelity * 100;
          case 4: return 100; // max out
          case 5: return layer.fibonacci_reliability * 100;
          default: return 0;
        }
      });

      const yScale = d3.scaleBand()
        .domain(layers.map((d) => d.name))
        .range([0, plotH])
        .padding(0.25);

      const xScale = d3.scaleLinear().domain([0, 100]).range([0, plotW]);

      // Grid
      g.append("g").selectAll("line")
        .data(xScale.ticks(5))
        .join("line")
        .attr("x1", (d) => xScale(d)).attr("x2", (d) => xScale(d))
        .attr("y1", 0).attr("y2", plotH)
        .attr("stroke", "#2a2a2a").attr("stroke-dasharray", "2,4");

      // Bars
      layers.forEach((layer, i) => {
        g.append("rect")
          .attr("x", 0)
          .attr("y", yScale(layer.name))
          .attr("width", xScale(barValues[i]))
          .attr("height", yScale.bandwidth())
          .attr("fill", layer.pass ? "#22c55e" : "#ef4444")
          .attr("opacity", 0.8)
          .attr("rx", 3);

        // PASS badge
        g.append("text")
          .attr("x", xScale(barValues[i]) + 8)
          .attr("y", yScale(layer.name) + yScale.bandwidth() / 2)
          .attr("dy", "0.35em")
          .attr("fill", layer.pass ? "#22c55e" : "#ef4444")
          .attr("font-size", "10px")
          .attr("font-weight", "700")
          .text(layer.pass ? "PASS" : "FAIL");

        // Metric label
        g.append("text")
          .attr("x", xScale(barValues[i]) + 42)
          .attr("y", yScale(layer.name) + yScale.bandwidth() / 2)
          .attr("dy", "0.35em")
          .attr("fill", "#ccc")
          .attr("font-size", "10px")
          .text(metricLabels[i]);
      });

      // Y axis (layer names)
      g.append("g")
        .call(d3.axisLeft(yScale).tickSize(0))
        .call((el) => el.select(".domain").remove())
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"));

      // Title
      g.append("text").attr("x", plotW / 2).attr("y", -18)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "14px")
        .attr("font-weight", "600").text("Full Computing Stack Validation");

      // Summary
      g.append("text").attr("x", plotW / 2).attr("y", plotH + 22)
        .attr("text-anchor", "middle").attr("fill", "#22c55e").attr("font-size", "11px")
        .text(`${data.overall.layers_passed}/${data.overall.total_layers} layers validated`);
    });
  }, [data]);

  return (
    <svg
      ref={svgRef}
      viewBox="0 0 800 300"
      className="w-full h-auto"
      style={{ background: "#1b1b1b", borderRadius: "8px" }}
    />
  );
};

export default ComputingStackChart;
