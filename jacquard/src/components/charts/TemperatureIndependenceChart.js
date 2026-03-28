import { useEffect, useRef, useState } from "react";

const TemperatureIndependenceChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetch("/data/temperature_velocity_independence.json")
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
      const margin = { top: 40, right: 20, bottom: 50, left: 55 };
      const midGap = 60;
      const leftWidth = (width - midGap) / 2;
      const rightWidth = (width - midGap) / 2;
      const plotH = height - margin.top - margin.bottom;

      // LEFT PANEL: Graph metrics vs temperature
      const leftG = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
      const leftPlotW = leftWidth - margin.left - 10;

      const temps = data.temperature_independence.temperatures;
      const metrics = data.temperature_independence.metrics;

      // Pick a few representative metrics
      const metricKeys = ["n_nodes", "n_edges", "degree_mean", "clustering_mean"];
      const metricColors = ["#58E6D9", "#B63E96", "#f59e0b", "#22c55e"];
      const metricLabels = ["Nodes", "Edges", "Degree Mean", "Clustering"];

      const xScale = d3.scaleLinear().domain([0, 11]).range([0, leftPlotW]);

      // Y axis - normalized 0-1
      const yScale = d3.scaleLinear().domain([0, 1.2]).range([plotH, 0]);

      // Grid
      leftG.append("g").selectAll("line")
        .data(yScale.ticks(4))
        .join("line")
        .attr("x1", 0).attr("x2", leftPlotW)
        .attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d))
        .attr("stroke", "#2a2a2a").attr("stroke-dasharray", "2,4");

      // X axis
      leftG.append("g")
        .attr("transform", `translate(0,${plotH})`)
        .call(d3.axisBottom(xScale).ticks(5))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      // Y axis
      leftG.append("g")
        .call(d3.axisLeft(yScale).ticks(4).tickFormat((d) => d.toFixed(1)))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      // Plot flat lines for each metric
      metricKeys.forEach((key, i) => {
        const vals = metrics[key].values;
        const maxVal = d3.max(vals) || 1;
        const normalized = vals.map((v) => v / maxVal);

        const line = d3.line()
          .x((d, j) => xScale(temps[j]))
          .y((d) => yScale(d));

        leftG.append("path")
          .datum(normalized)
          .attr("d", line)
          .attr("fill", "none")
          .attr("stroke", metricColors[i])
          .attr("stroke-width", 2);

        // Points
        normalized.forEach((v, j) => {
          leftG.append("circle")
            .attr("cx", xScale(temps[j]))
            .attr("cy", yScale(v))
            .attr("r", 3.5)
            .attr("fill", metricColors[i]);
        });
      });

      // Legend
      const legend = leftG.append("g").attr("transform", `translate(${leftPlotW - 110}, 0)`);
      metricKeys.forEach((key, i) => {
        legend.append("line")
          .attr("x1", 0).attr("x2", 14)
          .attr("y1", i * 16 + 4).attr("y2", i * 16 + 4)
          .attr("stroke", metricColors[i]).attr("stroke-width", 2);
        legend.append("text")
          .attr("x", 18).attr("y", i * 16 + 8)
          .attr("fill", "#ccc").attr("font-size", "9px").text(metricLabels[i]);
      });

      leftG.append("text").attr("x", leftPlotW / 2).attr("y", -18)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "13px")
        .attr("font-weight", "600").text("Temperature Independence");

      leftG.append("text").attr("x", leftPlotW / 2).attr("y", plotH + 38)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "10px")
        .text("Temperature (a.u.)");

      leftG.append("text").attr("x", leftPlotW / 2).attr("y", plotH + 26)
        .attr("text-anchor", "middle").attr("fill", "#22c55e").attr("font-size", "9px")
        .text("Variance = 0 across all temperatures");

      // RIGHT PANEL: Velocity blindness
      const rightG = svg.append("g").attr("transform", `translate(${leftWidth + midGap},${margin.top})`);
      const rightPlotW = rightWidth - 30;

      const pathsMatch = data.velocity_blindness.paths_match;
      const matchRate = data.velocity_blindness.match_rate;
      const nTrials = data.velocity_blindness.n_trials;

      // Show as a histogram of velocity differences, colored by match status
      const velDiffs = data.velocity_blindness.velocity_differences;
      const bins = d3.bin().domain([0, 8]).thresholds(12)(velDiffs);

      const xRight = d3.scaleLinear().domain([0, 8]).range([0, rightPlotW]);
      const yRight = d3.scaleLinear()
        .domain([0, d3.max(bins, (d) => d.length)])
        .range([plotH, 0]);

      // X axis
      rightG.append("g")
        .attr("transform", `translate(0,${plotH})`)
        .call(d3.axisBottom(xRight).ticks(4))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      // Y axis
      rightG.append("g")
        .call(d3.axisLeft(yRight).ticks(4))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      // Bars
      bins.forEach((bin) => {
        rightG.append("rect")
          .attr("x", xRight(bin.x0) + 1)
          .attr("y", yRight(bin.length))
          .attr("width", Math.max(0, xRight(bin.x1) - xRight(bin.x0) - 2))
          .attr("height", plotH - yRight(bin.length))
          .attr("fill", "#22c55e")
          .attr("opacity", 0.7)
          .attr("rx", 2);
      });

      rightG.append("text").attr("x", rightPlotW / 2).attr("y", -18)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "13px")
        .attr("font-weight", "600").text("Velocity Blindness");

      rightG.append("text").attr("x", rightPlotW / 2).attr("y", plotH + 26)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "10px")
        .text("Velocity difference");

      rightG.append("text").attr("x", rightPlotW / 2).attr("y", plotH + 38)
        .attr("text-anchor", "middle").attr("fill", "#22c55e").attr("font-size", "9px")
        .text(`${nTrials}/${nTrials} paths match (${(matchRate * 100).toFixed(0)}%)`);
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

export default TemperatureIndependenceChart;
