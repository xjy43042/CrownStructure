# -*- coding: utf-8 -*-
"""Resample 1m tree height tif files to 30m resolution NetCDF format."""
import os
from osgeo import gdal, osr
from multiprocessing import Pool
import subprocess

def raster_metadata(data):
    """Extract key metadata from raster file"""
    ds = gdal.Open(data)
    metadata = {
        'format': ds.GetDriver().ShortName,
        'band_count': ds.RasterCount,
        'nodata': ds.GetRasterBand(1).GetNoDataValue(),
        'band_stat': ds.GetRasterBand(1).GetStatistics(True, True),
        'proj4': osr.SpatialReference(ds.GetProjection()).ExportToProj4(),
        'wkt': ds.GetProjection(),
        'geocs': osr.SpatialReference(ds.GetProjection()).GetAttrValue('GEOGCS'),
        'uom': osr.SpatialReference(ds.GetProjection()).GetAttrValue('UNIT'),
        'projcs': osr.SpatialReference(ds.GetProjection()).GetAttrValue('PROJCS'),
        'epsg': osr.SpatialReference(ds.GetProjection()).GetAuthorityCode(None),
        'is_projected': osr.SpatialReference(ds.GetProjection()).IsProjected(),
        'geotransform': ds.GetGeoTransform(),
        'resolution': ds.GetGeoTransform()[1],
        'width': ds.RasterXSize * ds.GetGeoTransform()[1],
        'height': ds.RasterYSize * ds.GetGeoTransform()[5],
        'size_width': ds.RasterXSize,
        'size_height': ds.RasterYSize,
        'extent': [
            ds.GetGeoTransform()[0],
            ds.GetGeoTransform()[0] + ds.RasterXSize * ds.GetGeoTransform()[1],
            ds.GetGeoTransform()[3] + ds.RasterYSize * ds.GetGeoTransform()[5],
            ds.GetGeoTransform()[3]
        ],
        'centroid': [
            ds.GetGeoTransform()[0] + ds.RasterXSize * ds.GetGeoTransform()[1] / 2,
            ds.GetGeoTransform()[3] + ds.RasterYSize * ds.GetGeoTransform()[5] / 2
        ]
    }
    ds = None
    return metadata

def process_raster_file(filename, lltpath, lltpath2, llcpath):
    """
    Process single raster file through workflow:
    1. Set missing values
    2. Reproject to WGS84
    3. Convert to NetCDF
    """
    meta = raster_metadata(filename)
    s_srs = meta['proj4']
    t_srs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

    # Step 1: Set missing values
    command_setmiss = [
        'gdal_calc.py',
        '-A', filename,
        '--outfile', f'{lltpath}{filename}',
        '--calc', 'A*(A>=1)',
        '--NoDataValue', '0',  # missing values are 0
        '--overwrite'
    ]
    subprocess.run(command_setmiss, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Step 2: Reproject to geographic coordinates
    command_warp = [
        'gdalwarp',
        '-s_srs', s_srs,
        '-t_srs', t_srs,
        '-tr', '0.00025', '0.00025',
        '-r', 'max',
        '-srcnodata', '0',
        '-dstnodata', '0',
        '-co', 'COMPRESS=DEFLATE',
        f'{lltpath}{filename}',
        f'{lltpath2}{filename}'
    ]
    subprocess.run(command_warp, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Step 3: Convert to compressed NetCDF
    command_trans = [
        'gdal_translate',
        '-ot', 'Int16',
        '-of', 'netCDF',
        '-co', 'COMPRESS=DEFLATE',
        '-co', 'ZLEVEL=8',
        f'{lltpath2}{filename}',
        f'{llcpath}{os.path.splitext(filename)[0]}.nc'
    ]
    subprocess.run(command_trans, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Cleanup temporary files
    try:
        os.remove(lltpath + filename)
        os.remove(lltpath2 + filename)
    except OSError as e:
        print(f"Error: {e.strerror}")

def main():
    # Configure paths
    lltpath = '/tera10/yuanhua/xiangjy/CHM_Tolan/mkGlobe/new/20k/tif_setmiss1/'
    lltpath2 = '/tera10/yuanhua/xiangjy/CHM_Tolan/mkGlobe/new/20k/tif_tmp_30mmax/'
    llcpath = '/tera10/yuanhua/xiangjy/CHM_Tolan/mkGlobe/new/20k/nc_30mmax/'

    # Get input files
    allfile = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.tif')]

    # Process files in parallel
    max_processes = 20
    with Pool(processes=max_processes) as pool:
        pool.starmap(process_raster_file, [(filename, lltpath, lltpath2, llcpath) for filename in allfile])

if __name__ == "__main__":
    main()
