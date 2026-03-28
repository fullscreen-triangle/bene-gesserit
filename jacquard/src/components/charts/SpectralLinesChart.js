import { useEffect, useRef, useState } from "react";

const SpectralLinesChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetch("/data/hydrogen_spectral_lines.json")
      .then((r) => r.json())
      .then(setData);
  }, []);

  useEffect(() => {
    if (!data || !svgRef.current) return;

    import("d3").then((d3) => {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const width = 800;
      const height = 350;
      const margin = { top: 30, right: 30, bottom: 70, left: 60 };
      const topH = 200;
      const botH = 70;
      const gap = 20;
      const plotW = width - margin.left - margin.right;

      const results = data.results;
      const lines = results.map((d) => d.line.replace("_", " "));

      // Top panel: grouped bar chart
      const topG = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

      const x0 = d3.scaleBand().domain(lines).range([0, plotW]).padding(0.3);
      const x1 = d3.scaleBand().domain(["derived", "nist"]).range([0, x0.bandwidth()]).padding(0.08);
      const yMax = d3.max(results, (d) => Math.max(d.lambda_derived_nm, d.lambda_nist_nm));
      const y = d3.scaleLinear().domain([0, yMax * 1.1]).range([topH, 0]);

      // X axis
      topG.append("g")
        .attr("transform", `translate(0,${topH})`)
        .call(d3.axisBottom(x0).tickSize(0))
        .call((g) => g.select(".domain").attr("stroke", "#555"))
        .call((g) => g.selectAll(".tick text").attr("fill", "#ccc").attr("font-size", "10px")
          .attr("transform", "rotate(-20)").attr("text-anchor", "end"));

      // Y axis
      topG.append("g")
        .call(d3.axisLeft(y).ticks(5).tickFormat((d) => d + " nm"))
        .call((g) => g.select(".domain").attr("stroke", "#555"))
        .call((g) => g.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((g) => g.selectAll(".tick line").attr("stroke", "#444"));

      // Grid
      topG.append("g").selectAll("line")
        .data(y.ticks(5))
        .join("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", (d) => y(d)).attr("y2", (d) => y(d))
        .attr("stroke", "#333").attr("stroke-dasharray", "2,4");

      // Bars
      const colors = { derived: "#58E6D9", nist: "#B63E96" };

      results.forEach((d, i) => {
        const lineName = lines[i];
        const g = topG.append("g").attr("transform", `translate(${x0(lineName)},0)`);

        g.append("rect")
          .attr("x", x1("derived")).attr("y", y(d.lambda_derived_nm))
          .attr("width", x1.bandwidth()).attr("height", topH - y(d.lambda_derived_nm))
          .attr("fill", colors.derived).attr("opacity", 0.85).attr("rx", 2);

        g.append("rect")
          .attr("x", x1("nist")).attr("y", y(d.lambda_nist_nm))
          .attr("width", x1.bandwidth()).attr("height", topH - y(d.lambda_nist_nm))
          .attr("fill", colors.nist).attr("opacity", 0.85).attr("rx", 2);
      });

      // Legend
      const legend = topG.append("g").attr("transform", `translate(${plotW - 180}, 5)`);
      [{ label: "Derived", color: "#58E6D9" }, { label: "NIST", color: "#B63E96" }].forEach((item, i) => {
        legend.append("rect").attr("x", i * 90).attr("y", 0).attr("width", 12).attr("height", 12)
          .attr("fill", item.color).attr("rx", 2);
        legend.append("text").attr("x", i * 90 + 16).attr("y", 10)
          .attr("fill", "#e5e5e5").attr("font-size", "11px").text(item.label);
      });

      topG.append("text").attr("x", plotW / 2).attr("y", -12)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "13px")
        .attr("font-weight", "600").text("Hydrogen Spectral Lines: Derived vs NIST");

      // Bottom panel: error bars
      const botG = svg.append("g").attr("transform", `translate(${margin.left},${margin.top + topH + gap})`);

      const xBot = d3.scaleBand().domain(lines).range([0, plotW]).padding(0.3);
      const yBot = d3.scaleLinear().domain([0, 0.06]).range([botH, 0]);

      botG.append("g")
        .attr("transform", `translate(0,${botH})`)
        .call(d3.axisBottom(xBot).tickSize(0))
        .call((g) => g.select(".domain").attr("stroke", "#555"))
        .call((g) => g.selectAll(".tick text").attr("fill", "#ccc").attr("font-size", "9px")
          .attr("transform", "rotate(-20)").attr("text-anchor", "end"));

      botG.append("g")
        .call(d3.axisLeft(yBot).ticks(3).tickFormat((d) => (d * 100).toFixed(2) + "%"))
        .call((g) => g.select(".domain").attr("stroke", "#555"))
        .call((g) => g.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "9px"))
        .call((g) => g.selectAll(".tick line").attr("stroke", "#444"));

      results.forEach((d, i) => {
        const lineName = lines[i];
        botG.append("rect")
          .attr("x", xBot(lineName) + xBot.bandwidth() * 0.15)
          .attr("y", yBot(d.pct_error / 100))
          .attr("width", xBot.bandwidth() * 0.7)
          .attr("height", botH - yBot(d.pct_error / 100))
          .attr("fill", "#B63E96")
          .attr("opacity", 0.7)
          .attr("rx", 2);

        botG.append("text")
          .attr("x", xBot(lineName) + xBot.bandwidth() / 2)
          .attr("y", yBot(d.pct_error / 100) - 3)
          .attr("text-anchor", "middle")
          .attr("fill", "#e5e5e5")
          .attr("font-size", "8px")
          .text(d.pct_error.toFixed(3) + "%");
      });

      botG.append("text").attr("x", plotW / 2).attr("y", -6)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px")
        .text("Error (%)");
    });
  }, [data]);

  return (
    <svg
      ref={svgRef}
      viewBox="0 0 800 350"
      className="w-full h-auto"
      style={{ background: "#1b1b1b", borderRadius: "8px" }}
    />
  );
};

export default SpectralLinesChart;
