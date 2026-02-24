const path = require('path');
const CopyWebpackPlugin = require('copy-webpack-plugin');

module.exports = {
  entry: {
    'offscreen/offscreen-bundle': './src/offscreen.js',
  },
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: '[name].js',
    clean: true,
  },
  resolve: {
    extensions: ['.js'],
  },
  module: {
    rules: [
      {
        test: /\.js$/,
        exclude: /node_modules/,
        type: 'javascript/auto',
      },
    ],
  },
  plugins: [
    new CopyWebpackPlugin({
      patterns: [
        { from: 'manifest.json', to: '.' },
        { from: 'icons', to: 'icons' },
        { from: 'assets', to: 'assets' },
        { from: 'popup', to: 'popup' },
        { from: 'options', to: 'options' },
        { from: 'background', to: 'background' },
        { from: 'offscreen/offscreen.html', to: 'offscreen/offscreen.html' },
        { from: 'offscreen/posture-engine.js', to: 'offscreen/posture-engine.js' },
        { from: 'content', to: 'content' },
        { from: 'shared', to: 'shared' },
        { from: 'companion', to: 'companion' },
      ],
    }),
  ],
  // Don't split TF.js into chunks — extension needs single bundle
  optimization: {
    splitChunks: false,
  },
  performance: {
    maxAssetSize: 15 * 1024 * 1024, // TF.js is large
    maxEntrypointSize: 15 * 1024 * 1024,
  },
};
