import { useEffect, useRef, useState } from "react";

const CapacityChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);
  const [domainData, setDomainData] = useState(null);

  useEffect(() => {
    Promise.all([
      fetch("/data/shell_capacity.json").then((r) => r.json()),
      fetch("/data/cross_domain_capacity_validation.json").then((r) => r.json()),
    ]).then(([shell, domain]) => {
      setData(shell);
      setDomainData(domain);
    });
  }, []);

  useEffect(() => {
    if (!data || !domainData || !svgRef.current) return;

    import("d3").then((d3) => {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const width = 800;
      const height = 300;
      const margin = { top: 30, right: 20, bottom: 40, left: 50 };
      const midGap = 60;
      const leftWidth = (width - midGap) / 2;
      const rightWidth = (width - midGap) / 2;

      // Left panel: Line chart of C(n) = 2n^2
      const leftG = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
      const plotW = leftWidth - margin.left - 20;
      const plotH = height - margin.top - margin.bottom;

      const results = data.results;
      const xScale = d3.scaleLinear().domain([0.5, 7.5]).range([0, plotW]);
      const yScale = d3.scaleLinear().domain([0, 105]).range([plotH, 0]);

      // Grid lines
      leftG.append("g").selectAll("line")
        .data(yScale.ticks(5))
        .join("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d))
        .attr("stroke", "#333").attr("stroke-dasharray", "2,4");

      // X axis
      leftG.append("g")
        .attr("transform", `translate(0,${plotH})`)
        .call(d3.axisBottom(xScale).ticks(7).tickFormat(d3.format("d")))
        .call((g) => g.select(".domain").attr("stroke", "#555"))
        .call((g) => g.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "11px"))
        .call((g) => g.selectAll(".tick line").attr("stroke", "#555"));

      // Y axis
      leftG.append("g")
        .call(d3.axisLeft(yScale).ticks(5))
        .call((g) => g.select(".domain").attr("stroke", "#555"))
        .call((g) => g.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "11px"))
        .call((g) => g.selectAll(".tick line").attr("stroke", "#555"));

      // Continuous curve C(n) = 2n^2
      const curveData = d3.range(0.5, 7.6, 0.1).map((n) => ({ n, c: 2 * n * n }));
      const line = d3.line().x((d) => xScale(d.n)).y((d) => yScale(d.c)).curve(d3.curveMonotoneX);
      leftG.append("path")
        .datum(curveData)
        .attr("d", line)
        .attr("fill", "none")
        .attr("stroke", "#58E6D9")
        .attr("stroke-width", 2)
        .attr("opacity", 0.5);

      // Data points
      leftG.selectAll("circle")
        .data(results)
        .join("circle")
        .attr("cx", (d) => xScale(d.n))
        .attr("cy", (d) => yScale(d.C_derived))
        .attr("r", 5)
        .attr("fill", "#58E6D9")
        .attr("stroke", "#1b1b1b")
        .attr("stroke-width", 1.5);

      // Labels
      leftG.append("text").attr("x", plotW / 2).attr("y", -10)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "13px")
        .attr("font-weight", "600").text("C(n) = 2n\u00B2");

      leftG.append("text").attr("x", plotW / 2).attr("y", plotH + 35)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px").text("Shell number n");

      leftG.append("text").attr("transform", "rotate(-90)").attr("x", -plotH / 2).attr("y", -38)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px").text("Capacity C(n)");

      // Right panel: Horizontal bar chart
      const rightG = svg.append("g").attr("transform", `translate(${leftWidth + midGap},${margin.top})`);
      const rPlotW = rightWidth - 30;
      const domains = domainData.domains;
      const barH = plotH / domains.length;

      const xBarScale = d3.scaleLinear().domain([0, 100]).range([0, rPlotW]);
      const yBarScale = d3.scaleBand()
        .domain(domains.map((d) => d.domain))
        .range([0, plotH])
        .padding(0.25);

      // Bars
      rightG.selectAll("rect")
        .data(domains)
        .join("rect")
        .attr("x", 0)
        .attr("y", (d) => yBarScale(d.domain))
        .attr("width", (d) => xBarScale(d.pass_rate * 100))
        .attr("height", yBarScale.bandwidth())
        .attr("fill", "#58E6D9")
        .attr("opacity", 0.8)
        .attr("rx", 2);

      // Domain labels
      rightG.selectAll(".label")
        .data(domains)
        .join("text")
        .attr("class", "label")
        .attr("x", -5)
        .attr("y", (d) => yBarScale(d.domain) + yBarScale.bandwidth() / 2)
        .attr("dy", "0.35em")
        .attr("text-anchor", "end")
        .attr("fill", "#ccc")
        .attr("font-size", "9px")
        .text((d) => d.domain.replace(/\s*\(.*\)/, "").slice(0, 18));

      // Percentage labels
      rightG.selectAll(".pct")
        .data(domains)
        .join("text")
        .attr("class", "pct")
        .attr("x", (d) => xBarScale(d.pass_rate * 100) + 4)
        .attr("y", (d) => yBarScale(d.domain) + yBarScale.bandwidth() / 2)
        .attr("dy", "0.35em")
        .attr("fill", "#58E6D9")
        .attr("font-size", "10px")
        .attr("font-weight", "600")
        .text("100%");

      rightG.append("text").attr("x", rPlotW / 2).attr("y", -10)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "13px")
        .attr("font-weight", "600").text("Cross-Domain Validation");
    });
  }, [data, domainData]);

  return (
    <svg
      ref={svgRef}
      viewBox="0 0 800 300"
      className="w-full h-auto"
      style={{ background: "#1b1b1b", borderRadius: "8px" }}
    />
  );
};

export default CapacityChart;
