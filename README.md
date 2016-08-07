## pyCactus
You have reached the pyCactus documentation. This documentation is generated using sphinx. To view the documentation, download this branch and open `index.html` in a browser.

## Rebuilding documentation
1. Download and install sphinx.
2. Clone this repository 
```shell
git clone git@github.com:SNL-WaterPower/CACTUS-tools.git
```

3. Checkout the `doc` branch
```shell
git fetch
git checkout doc
```

4. Build the documentation locally.
```shell
make doc-build
```

5. If everything looks good and you'd like to update the documentation on the GitHub repository,
```shell
make doc-deploy
```


