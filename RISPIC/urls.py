"""RISPIC URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from . import views
from django.conf.urls import url
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    url(r'^$', views.index, name="index"),
    url(r'^signup/$', views.signup, name='signup'),
    url(r'^login', views.login, name='login'),
    url(r'^logout', views.logout, name='logout'),
    url(r'^pipeline', views.pipeline_start, name="pipeline"),
    url(r'^createresults', views.create_results, name="createresults"),
    url(r'^results', views.get_stored_results, name="results"),
    url(r'^readme', views.readme, name="readme"),
    url(r'^delete', views.delete, name="delete"),
    url(r'^remove', views.remove_result, name="remove"),
]
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)