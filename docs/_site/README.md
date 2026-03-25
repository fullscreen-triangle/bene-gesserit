# Bene Gesserit Documentation Site

This directory contains the source code for the **Bene Gesserit** documentation website, built with Jekyll and designed for GitHub Pages.

## 🌐 Live Site

Visit the live documentation at: [https://your-username.github.io/bene-gesserit](https://your-username.github.io/bene-gesserit)

## 🏗️ Site Structure

```
docs/
├── _config.yml              # Jekyll configuration
├── Gemfile                  # Ruby dependencies
├── index.md                 # Homepage
├── _fundamentals/           # Core concepts documentation
│   └── index.md
├── _language/               # Turbulance language documentation
│   └── index.md
├── _membrane_dynamics/      # Membrane computing documentation
│   └── index.md
├── _examples/               # Practical examples
│   └── index.md
└── README.md               # This file
```

## 🚀 Local Development

### Prerequisites

1. **Ruby** (version 2.7 or higher)
2. **Bundler** gem
3. **Git**

### Setup

1. **Install Ruby and Bundler** (if not already installed):
   ```bash
   # On macOS with Homebrew
   brew install ruby
   gem install bundler
   
   # On Ubuntu/Debian
   sudo apt-get install ruby-full build-essential zlib1g-dev
   gem install bundler
   
   # On Windows
   # Download and install Ruby from https://rubyinstaller.org/
   gem install bundler
   ```

2. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/bene-gesserit.git
   cd bene-gesserit/docs
   ```

3. **Install dependencies**:
   ```bash
   bundle install
   ```

4. **Serve the site locally**:
   ```bash
   bundle exec jekyll serve
   ```

5. **Open your browser** and navigate to `http://localhost:4000`

### Development Workflow

1. **Make changes** to any `.md` files or configuration
2. **Save the files** - Jekyll will automatically rebuild the site
3. **Refresh your browser** to see the changes
4. **Commit and push** when ready:
   ```bash
   git add .
   git commit -m "Update documentation"
   git push origin main
   ```

## 📝 Content Structure

### Homepage (`index.md`)
- Hero section with framework introduction
- Interactive Turbulance examples
- Quick start guide
- Feature highlights
- Navigation to main sections

### Fundamentals (`_fundamentals/`)
- Core concepts of biological quantum computing
- ATP-constrained dynamics
- Oscillatory entropy
- ENAQT mechanisms
- Mathematical foundations

### Turbulance Language (`_language/`)
- Complete language reference
- Syntax guide and examples
- Scientific method constructs
- Pattern matching and evidence collection
- Goal systems and metacognition

### Membrane Dynamics (`_membrane_dynamics/`)
- Biological quantum computing through membranes
- Maxwell's demons and information catalysis
- Circuit interface layers
- ENAQT transport mechanisms

### Examples (`_examples/`)
- Practical code examples
- Complete integration tutorials
- Scientific applications
- Step-by-step guides

## 🎨 Customization

### Styling
- CSS styles are embedded in each page for simplicity
- Modify the `<style>` sections in individual `.md` files
- Colors and themes can be customized in `_config.yml`

### Interactive Elements
- JavaScript for tabs and interactive demos
- Embedded in `<script>` sections of relevant pages
- Syntax highlighting via Rouge (configured in `_config.yml`)

### Navigation
- Main navigation defined in `_config.yml`
- Automatic breadcrumbs and section linking
- Responsive design for mobile devices

## 🔧 Configuration

### `_config.yml` Settings

```yaml
# Site information
title: "Bene Gesserit: Biological Quantum Computing Framework"
description: "Revolutionary framework combining ATP dynamics..."
url: "https://your-username.github.io"
baseurl: "/bene-gesserit"

# Build settings
markdown: kramdown
highlighter: rouge
theme: minima

# Collections for organized content
collections:
  fundamentals:
    output: true
    permalink: /:collection/:name/
  # ... other collections
```

### GitHub Pages Deployment

1. **Enable GitHub Pages** in your repository settings:
   - Go to Settings → Pages
   - Select "Deploy from a branch"
   - Choose "main" branch and "/docs" folder
   - Save

2. **Custom domain** (optional):
   - Add a `CNAME` file to the docs directory
   - Configure DNS settings with your domain provider

3. **Automatic deployment**:
   - GitHub Pages automatically rebuilds on every push to main
   - Check the Actions tab for build status

## 📱 Features

### Interactive Demos
- Tabbed interfaces for different examples
- Syntax-highlighted code blocks
- Responsive design for all devices

### Scientific Focus
- LaTeX math rendering (via MathJax if needed)
- Scientific notation and units
- Comprehensive code examples

### SEO Optimized
- Meta tags and descriptions
- Sitemap generation
- Social media cards

### Performance
- Optimized images and assets
- Minimal JavaScript for fast loading
- CDN-friendly static generation

## 🤝 Contributing

### Adding New Content

1. **Create new pages** in the appropriate collection directory
2. **Use the frontmatter format**:
   ```yaml
   ---
   layout: default
   title: "Page Title"
   permalink: /section/page-name/
   description: "Page description for SEO"
   ---
   ```

3. **Follow the established style**:
   - Use consistent heading structures
   - Include code examples
   - Add interactive elements where appropriate

### Reporting Issues

- Use GitHub Issues for bugs or content suggestions
- Include specific page URLs and descriptions
- Screenshots are helpful for layout issues

### Pull Requests

- Fork the repository
- Create a feature branch
- Make your changes
- Test locally with `bundle exec jekyll serve`
- Submit a pull request with clear description

## 📚 Resources

- [Jekyll Documentation](https://jekyllrb.com/docs/)
- [GitHub Pages Documentation](https://docs.github.com/en/pages)
- [Markdown Guide](https://www.markdownguide.org/)
- [Liquid Template Language](https://shopify.github.io/liquid/)

## 📄 License

This documentation is part of the Bene Gesserit project and follows the same license terms as the main repository.

---

**Happy documenting!** 🚀

For questions or support, please open an issue in the main repository or contact the development team. 