import { useEffect, useRef, useState } from "react";

const SpeedupChart = () => {
  const svgRef = useRef(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetch("/data/categorical_speedup_validation.json")
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

      // Combine measured + extrapolated
      const measured = data.per_size;
      const extrap = data.extrapolations;

      const allPoints = [
        ...measured.map((d) => ({
          N: d.N,
          cat: d.categorical_ops,
          conv: d.conventional_ops,
          speedup: d.speedup,
          extrapolated: false,
        })),
        ...extrap.map((d) => ({
          N: d.N_value,
          cat: d.categorical_ops_predicted,
          conv: d.conventional_ops_predicted,
          speedup: d.predicted_speedup,
          extrapolated: true,
        })),
      ];

      const xScale = d3.scaleLog().domain([50, 2e12]).range([0, plotW]);
      const yScale = d3.scaleLog().domain([1, 5e13]).range([plotH, 0]);

      // Grid
      g.append("g").selectAll("line")
        .data(yScale.ticks(6))
        .join("line")
        .attr("x1", 0).attr("x2", plotW)
        .attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d))
        .attr("stroke", "#2a2a2a").attr("stroke-dasharray", "2,4");

      g.append("g").selectAll("line")
        .data(xScale.ticks(6))
        .join("line")
        .attr("x1", (d) => xScale(d)).attr("x2", (d) => xScale(d))
        .attr("y1", 0).attr("y2", plotH)
        .attr("stroke", "#2a2a2a").attr("stroke-dasharray", "2,4");

      // Axes
      g.append("g")
        .attr("transform", `translate(0,${plotH})`)
        .call(d3.axisBottom(xScale).ticks(6, ".0e"))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      g.append("g")
        .call(d3.axisLeft(yScale).ticks(6, ".0e"))
        .call((el) => el.select(".domain").attr("stroke", "#555"))
        .call((el) => el.selectAll(".tick text").attr("fill", "#e5e5e5").attr("font-size", "10px"))
        .call((el) => el.selectAll(".tick line").attr("stroke", "#555"));

      // Conventional line
      const convLine = d3.line()
        .x((d) => xScale(d.N))
        .y((d) => yScale(d.conv))
        .curve(d3.curveMonotoneX);

      g.append("path")
        .datum(allPoints)
        .attr("d", convLine)
        .attr("fill", "none")
        .attr("stroke", "#B63E96")
        .attr("stroke-width", 2.5)
        .attr("stroke-dasharray", (d) => null);

      // Categorical line
      const catLine = d3.line()
        .x((d) => xScale(d.N))
        .y((d) => yScale(d.cat))
        .curve(d3.curveMonotoneX);

      g.append("path")
        .datum(allPoints)
        .attr("d", catLine)
        .attr("fill", "none")
        .attr("stroke", "#58E6D9")
        .attr("stroke-width", 2.5);

      // Measured points
      allPoints.filter((d) => !d.extrapolated).forEach((d) => {
        g.append("circle")
          .attr("cx", xScale(d.N)).attr("cy", yScale(d.conv))
          .attr("r", 5).attr("fill", "#B63E96").attr("stroke", "#1b1b1b").attr("stroke-width", 1.5);
        g.append("circle")
          .attr("cx", xScale(d.N)).attr("cy", yScale(d.cat))
          .attr("r", 5).attr("fill", "#58E6D9").attr("stroke", "#1b1b1b").attr("stroke-width", 1.5);
      });

      // Extrapolated points (hollow)
      allPoints.filter((d) => d.extrapolated).forEach((d) => {
        g.append("circle")
          .attr("cx", xScale(d.N)).attr("cy", yScale(d.conv))
          .attr("r", 5).attr("fill", "none").attr("stroke", "#B63E96").attr("stroke-width", 2);
        g.append("circle")
          .attr("cx", xScale(d.N)).attr("cy", yScale(d.cat))
          .attr("r", 5).attr("fill", "none").attr("stroke", "#58E6D9").attr("stroke-width", 2);
      });

      // Speedup annotations
      const annotate = [
        { N: 10000, speedup: "39x", yOff: -15 },
        { N: 1e6, speedup: "1.3M x", yOff: -15 },
        { N: 1e12, speedup: "1.3T x", yOff: -15 },
      ];
      annotate.forEach((a) => {
        const pt = allPoints.find((d) => d.N === a.N);
        if (!pt) return;
        const midY = (yScale(pt.cat) + yScale(pt.conv)) / 2;
        g.append("text")
          .attr("x", xScale(a.N) + 8)
          .attr("y", midY + a.yOff)
          .attr("fill", "#e5e5e5")
          .attr("font-size", "10px")
          .attr("font-weight", "600")
          .text(a.speedup);

        // Dashed connector
        g.append("line")
          .attr("x1", xScale(a.N)).attr("x2", xScale(a.N))
          .attr("y1", yScale(pt.cat)).attr("y2", yScale(pt.conv))
          .attr("stroke", "#666")
          .attr("stroke-dasharray", "3,3")
          .attr("stroke-width", 1);
      });

      // Legend
      const legend = g.append("g").attr("transform", `translate(${plotW - 220}, 10)`);
      [
        { label: "Conventional O(N log N)", color: "#B63E96" },
        { label: "Categorical O(log\u2083 N)", color: "#58E6D9" },
      ].forEach((item, i) => {
        legend.append("line")
          .attr("x1", 0).attr("x2", 20)
          .attr("y1", i * 20 + 5).attr("y2", i * 20 + 5)
          .attr("stroke", item.color).attr("stroke-width", 2.5);
        legend.append("circle")
          .attr("cx", 10).attr("cy", i * 20 + 5)
          .attr("r", 4).attr("fill", item.color);
        legend.append("text")
          .attr("x", 26).attr("y", i * 20 + 9)
          .attr("fill", "#e5e5e5").attr("font-size", "11px").text(item.label);
      });

      // Title
      g.append("text").attr("x", plotW / 2).attr("y", -18)
        .attr("text-anchor", "middle").attr("fill", "#e5e5e5").attr("font-size", "14px")
        .attr("font-weight", "600").text("Categorical Speedup: O(log\u2083 N) vs O(N log N)");

      // Axis labels
      g.append("text").attr("x", plotW / 2).attr("y", plotH + 42)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px").text("Problem size N");
      g.append("text").attr("transform", "rotate(-90)").attr("x", -plotH / 2).attr("y", -55)
        .attr("text-anchor", "middle").attr("fill", "#999").attr("font-size", "11px").text("Operations");
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

export default SpeedupChart;
