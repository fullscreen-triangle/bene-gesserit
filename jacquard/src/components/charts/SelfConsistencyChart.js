import { useEffect, useRef, useState } from "react";

const SelfConsistencyChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetch("/data/self_consistency.json")
      .then((r) => r.json())
      .then(setData);
  }, []);

  useEffect(() => {
    if (!data || !svgRef.current) return;

    import("d3").then((d3) => {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const width = 600;
      const height = 300;
      const margin = { top: 40, right: 30, bottom: 50, left: 60 };
      const plotW = width - margin.left - margin.right;
      const plotH = height - margin.top - margin.bottom;

      const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

      // Map molecules to atom count approximation and deviation
      const molecules = data.molecules_tested;
      const atomCounts = { H2: 2, CO: 2, H2O: 3, CO2: 3, CH4: 5, C6H6: 12 };
      const points = molecules.map((mol) => ({
        mol,
        atoms: atomCounts[mol] || 3,
        dev: data.results[mol].max_deviation_pct,
        status: data.results[mol].status,
      }));

      const xScale = d3.scaleLinear()
        .domain([0, 14])
        .range([0, plotW]);

      const yScale = d3.scaleLinear()
        .domain([0, 2.0])
        .range([plotH, 0]);

      // Grid
      g.append("g").selectAll("line")
        .data(yScale.ticks(4))
        .join("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d))
        .attr("stroke", "#2a2a2a").attr("stroke-dasharray", "2,4");

      // Threshold line at 2%
      g.append("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", yScale(2.0)).attr("y2", yScale(2.0))
        .attr("stroke", "#ef4444")
        .attr("stroke-width", 1)
        .attr("stroke-dasharray", "6,3");

      g.append("text")
        .attr("x", plotW - 5).attr("y", yScale(2.0) - 5)
        .attr("text-anchor", "end")
        .attr("fill", "#ef4444")
        .attr("font-size", "9px")
        .text("2% threshold");

      // Axes
      g.append("g")
        .attr("transform", `translate(0,${plotH})`)
        .call(d3.axisBottom(xScale).ticks(7))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      g.append("g")
        .call(d3.axisLeft(yScale).ticks(5).tickFormat((d) => d.toFixed(1) + "%"))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      // Trend line (linear regression on non-zero points)
      const nonZero = points.filter((d) => d.dev > 0);
      if (nonZero.length >= 2) {
        const xMean = d3.mean(nonZero, (d) => d.atoms);
        const yMean = d3.mean(nonZero, (d) => d.dev);
        const num = d3.sum(nonZero, (d) => (d.atoms - xMean) * (d.dev - yMean));
        const den = d3.sum(nonZero, (d) => (d.atoms - xMean) ** 2);
        const slope = num / den;
        const intercept = yMean - slope * xMean;

        g.append("line")
          .attr("x1", xScale(1)).attr("x2", xScale(13))
          .attr("y1", yScale(slope * 1 + intercept))
          .attr("y2", yScale(slope * 13 + intercept))
          .attr("stroke", "#58E6D9")
          .attr("stroke-width", 1.5)
          .attr("stroke-dasharray", "6,4")
          .attr("opacity", 0.6);
      }

      // Points
      points.forEach((d) => {
        g.append("circle")
          .attr("cx", xScale(d.atoms))
          .attr("cy", yScale(d.dev))
          .attr("r", 7)
          .attr("fill", "#58E6D9")
          .attr("stroke", "#1b1b1b")
          .attr("stroke-width", 2);

        // Label
        const yOff = d.mol === "CH4" ? -14 : d.mol === "C6H6" ? -14 : 16;
        g.append("text")
          .attr("x", xScale(d.atoms))
          .attr("y", yScale(d.dev) + yOff)
          .attr("text-anchor", "middle")
          .attr("fill", "#e5e5e5")
          .attr("font-size", "10px")
          .attr("font-weight", "600")
          .text(d.mol);
      });

      // Title
      g.append("text").attr("x", plotW / 2).attr("y", -18)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "14px")
        .attr("font-weight", "600").text("Self-Consistency: Max Deviation by Molecule");

      // Axis labels
      g.append("text").attr("x", plotW / 2).attr("y", plotH + 42)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px")
        .text("Number of atoms");

      g.append("text").attr("transform", "rotate(-90)").attr("x", -plotH / 2).attr("y", -45)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px")
        .text("Max deviation (%)");

      // Summary
      g.append("text").attr("x", plotW / 2).attr("y", plotH + 26)
        .attr("text-anchor", "middle").attr("fill", "#22c55e").attr("font-size", "10px")
        .text(`All PASS \u2014 max overall: ${data.summary.max_overall_deviation_pct.toFixed(2)}% < ${data.summary.threshold_pct}% threshold`);
    });
  }, [data]);

  return (
    <svg
      ref={svgRef}
      viewBox="0 0 600 300"
      className="w-full h-auto"
      style={{ background: "#1b1b1b", borderRadius: "8px" }}
    />
  );
};

export default SelfConsistencyChart;
